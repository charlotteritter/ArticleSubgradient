$ONTEXT
The full "true" problem.
Check which scenario file are you importing

Production limits adapted via Morales  Latorre Ramos (2013)
x here is the production ABOVE the min g, so total production is g*r + x

$OFFTEXT

OPTIONS PROFILE =3, RESLIM   = 4000, LIMROW   = 5, LP = CPLEX, MIP = gurobi, RMIP=Gurobi, NLP = CONOPT, MINLP = DICOPT, MIQCP = CPLEX,
         SOLPRINT = OFF, decimals = 8, optcr=0.01, optca=0.01, threads =8, integer4=0;

********************************************************************************
* run_time_total recorder
SCALAR start_time, end_time, run_time_total;

*-------------------------------------------------------------------------------

** sets defined in input file
SETS T times/t1*t24/;
SETS W scenarios /scen1*scen1500/;

ALIAS (T,TT);
ALIAS (W,I);
** define generator costs and wind selling prices
TABLE COSTS(T,*)
$ONDELIM
$INCLUDE wind_costs.csv
$OFFDELIM
;

** define wind realizations at all time periods
TABLE WIND(W,T)
$ONDELIM
$INCLUDE wind_scenarios_150.csv
$OFFDELIM
;


scalar PROBABILITY;
PROBABILITY = 1/CARD(W);
;


** define tolerance threshold
SCALAR TOL;
TOL =0.05;

* find Ntol + 1st value
parameter maxwind(t), minwind(t), dummywind(w,t) ;
maxwind(t) =smax(w,wind(w,t)) ;
dummywind(w,t) = wind(w,t) ;

scalar it ;
it = floor(card(w)*tol) + 1;

* index of it
set dummy(w);
* make the dum_iter go till at least the size of it
set dum_iter /dum_iter1*dum_iter100/;
loop(t,
loop(dum_iter$(ord(dum_iter)le it),
* find the smallest wind value for this t
         minwind(t) = smin(w,dummywind(w,t)) ;
* index of smallest wind value
         dummy(w) = yes$(dummywind(w,t) eq minwind(t)) ;
* make the smallest wind value large
         dummywind(w,t)$dummy(w) =maxwind(t) ;
); );
display minwind ;

** define generator ramp constraints \Delta
* this data is from Ostrowski's paper, our generator was too too large
SCALAR ramp;
ramp=50;

** define upper bound of generator capacity
SCALAR G, GG, start_cost;
** Max Generator
G=130;
** Min Generator
GG=20; 

parameter BigMM(w), BigM(w,t);
BigMM(w)= smax(t, wind(w,t));
BigM(w,t)= G - wind(w,t) + minwind(t);



SCALARS UPTIME, DOWNTIME;
UPTIME =3;
DOWNTIME=3;

********************************************************************************
*                                begin model
********************************************************************************
POSITIVE VARIABLES X(W,T), Y(T) ;
VARIABLES OBJ;
BINARY VARIABLE Z(W), R(W,T),U(W,T), V(W,T) ;

EQUATIONS
        Objective
        Const1_1(W,T) bigm constraint
        Const1_2     sum probabilities
        Const_3_1(W,T) ramping constraint
        Const_3_2(W,T) ramping constraint
        Const_4_1(W,T) generator off constraints
        Const_4_2(W,T) generator off constraints

        Const_5(W,T) generator running constraints

        Const_6(W,T)

        ;

Objective.. OBJ=E= SUM(T,COSTS(T, 'REW')*Y(T) - PROBABILITY*SUM(W, COSTS(T, 'GEN')*( X(W,T) + GG*R(w,t) ) )  )  ;

*Const1_1(W,T).. Y(T)-X(W,T)-R(W,T)*GG -Z(W)*WIND(W,T) =L= (1-Z(W))*BigMM(w) ;
Const1_1(W,T).. Y(T)-X(W,T)-R(W,T)*GG -WIND(W,T) =L= Z(W)*BigM(w,t) ;

Const1_2..  SUM(W,Z(W)) =L= floor(card(w)*TOL) ;

* The generator constraints

*Ramp

Const_3_1(W,T)$( ord(t) le (card(t)-1)).. X(W,T+1) + GG*R(w,t+1)- X(W,T) - GG*R(w,t)  =L= ramp*(U(w,t+1) + R(w,t)) ;
Const_3_2(W,T)$( ord(t) le (card(t)-1)).. X(W,T) + GG*R(w,t)- X(W,T+1) -GG*R(w,t+1)=L= ramp*(V(w,t+1) + R(w,t+1)) ;

*On/off

Const_4_1(W,T)$(ord(t) ge uptime)..   Sum(tt$((ORD(TT) le ORD(T)) and (ORD(TT) ge (ORD(T) - UPTIME +1))),U(W,TT)) =L= R(W,T) ;
Const_4_2(W,T)$(ord(t) ge downtime )..Sum(tt$((ORD(TT) le ORD(T)) and (ORD(TT) ge (ORD(T) - DOWNTIME +1))),V(W,TT)) =L= 1-R(W,T) ;


Const_5(W,T)$(ord(t) ge 2).. U(W,T) - V(W,T) =E= R(W,T) - R(W,T-1) ;

*if uptime more than 1
Const_6(W,T)$( ord(t) le (card(t)-1))..   (G-GG)*R(W,T)- (G-ramp)*U(W,T) - (G-ramp)*V(W,T+1) =G= X(W,T);


*** bounds on any variables
X.UP(W,T) =  G-GG;
U.UP(W,T) = 1;
V.UP(W,T) =1;
Z.UP(w)=1;

* initialize the on/off variables
* assume the unit was on for last (uptime -1) periods
r.fx(w,t)$(ord(t) eq 1) =1;
u.fx(w,t)$(ord(t) eq 1) = 0;
v.fx(w,t)$(ord(t) eq 1) = 0 ;
* assume the generator was producing minimim power in last time period
x.up(w,t)$(ord(t) eq 1)= ramp - GG;



option limrow = 10000 ;
MODEL  SCHEDULE    /ALL/ ;
schedule.optfile =1;

start_time=jnow;
SOLVE SCHEDULE USING MIP MAXIMIZING OBJ;
display obj.l ;

end_time = jnow;
run_time_total = ghour(end_time - start_time)*3600 + gminute(end_time - start_time)*60 + gsecond(end_time - start_time);
display run_time_total, u.l, r.l, v.l,x.l, y.l ;


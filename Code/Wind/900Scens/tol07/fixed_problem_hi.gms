$ONTEXT
This is Step 5 of Algorithm 1 of file v10.pdf (iEEE paper)
The 1500 scenario fixed problem
Follow up of SAA.gms

$OFFTEXT

OPTIONS PROFILE =3, RESLIM   = 1800, LIMROW   = 5, LP = CPLEX, MIP = gurobi, RMIP=Gurobi, NLP = CONOPT, MINLP = DICOPT, MIQCP = CPLEX,
         SOLPRINT = OFF, decimals = 8, optcr=0.0, optca=0.0, threads =8, integer4=1;

********************************************************************************
* run_time_total recorder
SCALAR start_time, end_time, run_time_total;

*-------------------------------------------------------------------------------

** sets defined in input file
SETS T times/t1*t24/;
SETS W scenarios /scen1*scen900/;

ALIAS (T,TT);
ALIAS (W,I);
ALIAS (W,SCEN);
** define generator costs and wind selling prices
TABLE COSTS(T,*)
$ONDELIM
$INCLUDE wind_costs.csv
$OFFDELIM
;

** define wind realizations at all time periods
TABLE WIND(W,T)
$ONDELIM
$INCLUDE wind_scenarios.csv
$OFFDELIM
;

scalar PROBABILITY;
PROBABILITY = 1/CARD(W);
;


** define tolerance threshold
SCALAR TOL;
TOL =0.07;


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

** define generator ramp constraints
* this data is from Ostrowski's paper, our generator was too too large
SCALAR ramp;
ramp=50;

** define upper bound of generator capacity
SCALAR G, GG, start_cost;
G=130;
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
POSITIVE VARIABLES X(W,T), Y(T),U(W,T), V(W,T) ;
VARIABLES OBJ;
BINARY VARIABLE Z(W), R(W,T) ;

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

Const1_1(W,T).. Y(T)-X(W,T)-R(W,T)*GG -WIND(W,T) =L= Z(W)*BigM(w,t) ;

Const1_2..  SUM(W,Z(W)) =L= card(w)*TOL ;

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

set iter /iter1*iter30/;
parameter  profit(iter), y_previous(t), run_time(iter);
scalar profit_orig, t1, t2;


TABLE y_100(t,iter)
$ONDELIM
$INCLUDE sampled_dynamic_hi.csv
$OFFDELIM
;

MODEL  SCHEDULE    /ALL/ ;

* 1st solve with the plain old 20 scenario promise
* use a solution from 20 scenarios scenarios
*y.fx(t)=y_20(t,'values');
*SOLVE SCHEDULE USING MIP MAXIMIZING OBJ;
*profit_orig= obj.l;

*2nd iteration onward
loop(iter,
         y.fx(t)=y_100(t,iter);
         t1=jnow ;
         SOLVE SCHEDULE USING MIP MAXIMIZING OBJ;
         t2=jnow;
         run_time(iter) = ghour(t2 - t1)*3600 + gminute(t2 - t1)*60 + gsecond(t2 - t1);
         profit(iter)= obj.l;
);

display y.l;

********************************************************************************
*                                write output
********************************************************************************


FILE fixed_profit /fixed_profit_hi_2.csv/;
fixed_profit.PC = 5;
fixed_profit.ND = 3;
PUT fixed_profit;
loop(iter, put iter.tl put profit(iter) put run_time(iter) put /; );
PUTCLOSE fixed_profit;

scalar sma;
sma=smax(iter,profit(iter));
File TestingFile3 / Alg_hi.csv /;
TestingFile3.pc=5;
TestingFile3.nd=5;
put TestingFile3; 
put 'Omega', put 'Tolerance', put 'Step Size Rule', put 'Iterations', put 'Converged?', put 'Gap LR', put 'Gap Naive', put 'Obj. Naive', put 'Obj. LR', put 'Gap' put 'Time Naive', put 'Time LR', put 'Final Lambda', put 'LB Heuristic' put /;
put '', put '', put '', put '', put '', put '', put '', put '', put '', put '', put '', put '', put '', put sma put /;


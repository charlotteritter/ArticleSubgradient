$ONTEXT
This is the GAMS code corresponding to the paper "Statistical performance of subgradient step-size update rules in Lagrangian relaxations of chance-constrained optimization models" by C. Ritter and B. Singh
See: https://github.com/charlotteritter/ArticleSubgradient

This is Step 4 and 6 of Algorithm 2 for wind model (Model II)
After this file, we run fixed scenario problem - fixed_problem.gms to finish Algorithm 2
$OFFTEXT

OPTIONS PROFILE =3, RESLIM   = 2100, LIMROW   = 5, LP = CPLEX, MIP = gurobi, RMIP=Gurobi, NLP = CONOPT, MINLP = DICOPT, MIQCP = CPLEX,
         SOLPRINT = OFF, decimals = 8, optcr=0.0005, optca=0.0005, threads =8, integer4=1;

********************************************************************************
* run_time_total recorder
SCALAR start_time, end_time, run_time_total;

*-------------------------------------------------------------------------------

SETS T times/t1*t24/;
SETS W scenarios /w1*w7200/;
SETS iter iterations /iter1*iter30/;

SETS rho_ind /r1/;
parameter rho_val(rho_ind) / r1 0.1/;

ALIAS (T,TT);
ALIAS (W,I);
** define generator costs and wind selling prices
TABLE COSTS(T,*)
$ONDELIM
$INCLUDE wind_costs.csv
$OFFDELIM
;

** define wind realizations at all time periods.
TABLE WIND(W,T)
$ONDELIM
$INCLUDE wind_scenarios_7200.csv
$OFFDELIM
;

** all scenarios are weighed equally here
table probability(w,*)
$ONDELIM
$INCLUDE ScenarioWeights_naive.csv
$OFFDELIM
;

** define tolerance threshold
SCALAR TOL;
TOL =0.03;

* find Ntol + 1st value
parameter maxwind(t), minwind(t), dummywind(w,t) ;
scalar it ;

* index of it
set dummy_set(w);
* make the dum_iter go till at least the size of it
set dum_iter /dum_iter1*dum_iter100/;

parameter BigM(w,t);

** define generator ramp constraints
* this data is from Ostrowski's paper, our generator was too too large
SCALAR ramp;
ramp=50;

** define upper bound of generator capacity
SCALAR G, GG, start_cost;
G=130;
GG=20;


SCALARS UPTIME, DOWNTIME;
UPTIME =3;
DOWNTIME=3;

scalars m,n,size, step, delta, rho, new_n, new_m, tot_time ;    ;
parameters time(rho_ind,iter), promise(rho_ind,iter,t), profit(rho_ind,iter) , num_scen(rho_ind,iter), y_previous(t) ;

* number of initial scenarios and % increase at each iterations
step = 20;
delta = 0.1 ;




********************************************************************************
*                                begin model
********************************************************************************
POSITIVE VARIABLES X(W,T), Y(T),U(W,T), V(W,T), dummy(t) ;
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

        Const_dum1(t)
        Const_dum2(t)

        ;

Objective.. OBJ=E= SUM(T,COSTS(T, 'REW')*Y(T) - SUM(w$((ord(w) ge n) and (ord(w) le m)), probability(w,'value') *COSTS(T, 'GEN')*( X(W,T) + GG*R(w,t) ) )  )  - rho*sum(t,dummy(t)) ;

Const1_1(W,T)$((ord(w) ge n) and (ord(w) le m)).. Y(T)-X(W,T)-R(W,T)*GG -WIND(W,T) =L= Z(W)*BigM(w,t) ;

Const1_2..  SUM(W$((ord(w) ge n) and (ord(w) le m)),Z(W)*probability(w,'value')) =L= TOL ;


* The generator constraints

*Ramp

Const_3_1(W,T)$( (ord(t) le (card(t)-1)) and (ord(w) ge n) and (ord(w) le m)).. X(W,T+1) + GG*R(w,t+1)- X(W,T) - GG*R(w,t)  =L= ramp*(U(w,t+1) + R(w,t)) ;
Const_3_2(W,T)$( (ord(t) le (card(t)-1)) and (ord(w) ge n) and (ord(w) le m)).. X(W,T) + GG*R(w,t)- X(W,T+1) -GG*R(w,t+1)=L= ramp*(V(w,t+1) + R(w,t+1)) ;

*On/off

Const_4_1(W,T)$( (ord(t) ge uptime)   and (ord(w) ge n) and (ord(w) le m))..   Sum(tt$((ORD(TT) le ORD(T)) and (ORD(TT) ge (ORD(T) - UPTIME +1))),U(W,TT)) =L= R(W,T) ;
Const_4_2(W,T)$( (ord(t) ge downtime) and (ord(w) ge n) and (ord(w) le m))..   Sum(tt$((ORD(TT) le ORD(T)) and (ORD(TT) ge (ORD(T) - DOWNTIME +1))),V(W,TT)) =L= 1-R(W,T) ;


Const_5(W,T)$(  (ord(t) ge 2) and (ord(w) ge n) and (ord(w) le m)).. U(W,T) - V(W,T) =E= R(W,T) - R(W,T-1) ;

*if uptime more than 1
Const_6(W,T)$(  (ord(t) le (card(t)-1)) and (ord(w) ge n) and (ord(w) le m))..   (G-GG)*R(W,T)- (G-ramp)*U(W,T) - (G-ramp)*V(W,T+1) =G= X(W,T);

Const_dum1(t).. dummy(t)=g= y(t) - y_previous(t);
Const_dum2(t).. dummy(t)=g= y_previous(t) - y(t);



*** bounds on any variables
X.UP(W,T) =  G-GG;
U.UP(W,T) = 1;
V.UP(W,T) = 1;
Z.UP(w)   = 1;


MODEL  SCHEDULE    /ALL/ ;



********************************************************************************
*                                begin iterations
********************************************************************************

loop(rho_ind,   
tot_time = 0;
size =0;
*We take the first 20 scenarios 
new_n = 1;
new_m = 20;
y_previous(t) = 0;

loop(iter,

if (tot_time < 1800,

         rho =rho_val(rho_ind);
* at the first iteration no proximal term due to first solving the model (1) without the regularization
if (ord(iter) eq 1,  rho =0 ; options optca =0.04, optcr =0.04) ;
if (ord(iter) gt 1,  options optca =0.0005, optcr =0.0005) ;

* choose scenarios between n and m
         n= new_n   ;
         m= new_m ;
         size = m-n+1 ;



*****************************************************************
* choosing Big M
*****************************************************************

maxwind(t) =smax(w$((ord(w) ge n) and (ord(w) le m)),wind(w,t)) ;
dummywind(w,t) = wind(w,t) ;

it = floor(size*tol) + 1;


loop(t,
loop(dum_iter$(ord(dum_iter)le it),
* find the smallest wind value for this t
         minwind(t) = smin(w$((ord(w) ge n) and (ord(w) le m)),dummywind(w,t)) ;
* index of smallest wind value
         dummy_set(w)$((ord(w) ge n) and (ord(w) le m)) = yes$(dummywind(w,t) eq minwind(t)) ;
* make the smallest wind value large
         dummywind(w,t)$(dummy_set(w) and ((ord(w) ge n) and (ord(w) le m))) =maxwind(t) ;
); );

BigM(w,t)$((ord(w) ge n) and (ord(w) le m))= G - wind(w,t) + minwind(t);


******************************************************************
*start of bounds
******************************************************************

* initialize the on/off variables
* assume the unit was on for last (uptime -1) periods
r.fx(w,t)$((ord(t) eq 1) and (ord(w) ge n) and (ord(w) le m)) =1;
u.fx(w,t)$((ord(t) eq 1) and (ord(w) ge n) and (ord(w) le m)) = 0;
v.fx(w,t)$((ord(t) eq 1) and (ord(w) ge n) and (ord(w) le m)) = 0 ;
* assume the generator was producing minimim power in last time period
x.up(w,t)$((ord(t) eq 1) and (ord(w) ge n) and (ord(w) le m)) = ramp - GG;


****************************** end  of bounds

         start_time=jnow;
         SOLVE SCHEDULE USING MIP MAXIMIZING OBJ ;
         end_time = jnow;
         time(rho_ind,iter) =  ghour(end_time - start_time)*3600 + gminute(end_time - start_time)*60 + gsecond(end_time - start_time);
         promise(rho_ind,iter,t) = y.l(t) ;
         y_previous(t) = y.l(t) ;
         profit(rho_ind,iter) = obj.l ;
         num_scen(rho_ind,iter) = size ;

         new_n = m + 1;
         new_m= m + ceil(size*(1+delta)) ;
);

* calculate total time loop has been running
         tot_time = time(rho_ind,iter) + tot_time ;

);

);

display tot_time ;
display time, num_scen, profit;

FILE sampled_dynamic /sampled_dynamic.csv/;
sampled_dynamic.PC = 5;
sampled_dynamic.ND = 3;
PUT sampled_dynamic;
loop(rho_ind,LOOP(t, loop(iter, put promise(rho_ind,iter,t);) put /; ); );

PUTCLOSE sampled_dynamic;


FILE profit_value /profit_value.csv/;
profit_value.PC = 5;
profit_value.ND = 3;
PUT profit_value;
loop(rho_ind, put rho_val(rho_ind) put /;
put ''  put 'time' put ' number scenarios' put 'profit_value' put /;
LOOP(iter, put iter.tl PUT time(rho_ind,iter) put num_scen(rho_ind,iter) put profit(rho_ind,iter) put /;);
);
PUTCLOSE profit_value;



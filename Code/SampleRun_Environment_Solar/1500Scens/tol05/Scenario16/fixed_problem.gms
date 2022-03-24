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
SETS W scenarios /scen1*scen1500/;

ALIAS (T,TT);
ALIAS (W,I);
ALIAS (W,SCEN);

SETS iter iterations /iter1*iter30/;

TABLE Solar(scen,t)
$ondelim
*$INCLUDE %SOLAR%.csv
$INCLUDE solar_scenarios.csv
$offdelim
;

*Tolerance 
scalar tol;
*tol=%TOL%;
tol=0.05;


* time limit for each problem
scalar time_limit;
*time_limit=%TIMELIM%;
time_limit=2250;

ALIAS (T,TT);
alias(scen,i);

scalar n;
n=card(scen);

TABLE y_100(t,iter)
$ONDELIM
$INCLUDE sampled_dynamic.csv
$OFFDELIM
;


** define battery  operation costs costs and solar selling prices

TABLE PRICES(t,*)
$ONDELIM
$INCLUDE battery_revenue.csv
$OFFDELIM
;

Prices(t,'rew')     =  - Prices(t,'rew');
Prices(t,'char')    =  - Prices(t,'char');
Prices(t,'dischar') =  - Prices(t,'dischar');
** define solar scenarios at all time periods


* Scaling of Solar power scenarios ;
scalar scale ;
scale = 1;
Solar(scen,t) = scale* Solar(scen,t) ;
* Remove too many decimals in Solar
Solar(scen,t) = round(Solar(scen,t),2) ;


scalar PROBABILITY;
PROBABILITY = 1/CARD(W);
;

scalar eta ;
*from Ben paper
eta = 0.9
;
parameters max_store(t), min_store(t), max_charge, max_discharge;


** define tolerance threshold
SCALAR threshold;
threshold = floor(card(scen)*TOL)  ;

parameter BigX, LowX, X_0 maximum minimum initial energy stored ;
parameter BigM(scen,t) find a good BigM ;

BigX = 960 ;
LowX = 0.2* BigX ;
X_0  = 0.5* BigX ;
max_charge =  0.5* BigX ;
max_discharge =  0.5* BigX ;

************** Find a Big M
* find Ntol + 1st value
parameter maxsolar(t), minsolar(t), dummysolar(scen,t) ;
maxsolar(t) =smax(scen,solar(scen,t)) ;
dummysolar(scen,t) = solar(scen,t) ;

scalar it ;
it = floor(card(scen)*tol) + 1;

* index of it
set dummy_set(scen);
* make the dum_iter go till at least the size of it
set dum_iter /dum_iter1*dum_iter100/;
loop(t,
loop(dum_iter$(ord(dum_iter)le it),
* find the smallest solar value for this t
         minsolar(t) = smin(scen,dummysolar(scen,t)) ;
* index of smallest solar value
         dummy_set(scen) = yes$(dummysolar(scen,t) eq minsolar(t)) ;
* make the smallest solar value large
         dummysolar(scen,t)$dummy_set(scen) =maxsolar(t) ;
); );
scalar G upper bound on q - p ;
G = min(eta*(BigX - LowX), max_discharge) ;

BigM(scen,t)= G - solar(scen,t) + minsolar(t);

********************************************************************************
*                                begin model
********************************************************************************



POSITIVE VARIABLES P(scen,t), Q(scen,t), Y(T), X(scen,t) ;
VARIABLES OBJ;
BINARY VARIABLE Z(scen) ;

scalar counter ;

EQUATIONS
        Objective
        Const1(scen,t)    balance constraint
        Const2(scen,t)    max charge
        Const3(scen,t)    max discharge
        Const_chance_1(scen,t)    chance constraint big M
        Const_chance_PH(scen,t)
        Const_chance_2            chance constraint sum probabilities
        ;

Objective.. OBJ=E= SUM(T,Prices(T, 'REW')*Y(T) - PROBABILITY*Sum(w, ( Prices(T, 'CHAR')* P(w,t) + Prices(t, 'DISCHAR') * Q(w,t) ) ) )     ;

Const1(scen,t)$(ord(t) lt card(t))..
         X(scen,t+1) =E= X(scen,t) + eta* P(scen,t) - (1/eta)* Q(scen,t) ;

Const_chance_1(scen,t).. Y(T) + P(scen,t) -  Q(scen,t) -SOLAR(scen,t) =L= Z(scen)*BigM(scen,t) ;


Const_chance_2..      - sum(scen, z(scen)) =G= -threshold;


*** bounds on any variables
x.up(scen,t) = BigX ;
x.lo(scen,t) = LowX ;
q.up(scen,t) = max_discharge ;
p.up(scen,t) =  max_charge ;
x.fx(scen,'t1') = X_0 ;
z.prior(scen)   = 1;


parameter last_x(scen,t), last_p(scen,t), last_q(scen,t), last_z(scen), last_ph(scen) ;
******* ALL MODELS

model schedule     / Objective,  Const1, Const_chance_1, Const_chance_2/ ;

parameter  profit(iter), y_previous(t), run_time(iter);
scalar profit_orig, t1, t2,start_time,end_time,tot_time;

start_time=jnow;
loop(iter,
         y.fx(t)=y_100(t,iter);
         t1=jnow ;
         SOLVE SCHEDULE USING MIP MINIMIZING OBJ;
         t2=jnow;
         run_time(iter) = ghour(t2 - t1)*3600 + gminute(t2 - t1)*60 + gsecond(t2 - t1);
         profit(iter)= obj.l;
);
end_time=jnow;
display y.l;

tot_time =  ghour(end_time - start_time)*3600 + gminute(end_time - start_time)*60 + gsecond(end_time - start_time);


********************************************************************************
*                                write output
********************************************************************************


FILE fixed_profit /fixed_profit.csv/;
fixed_profit.PC = 5;
fixed_profit.ND = 3;
PUT fixed_profit;
loop(iter, put iter.tl put profit(iter) put run_time(iter) put /; );
PUTCLOSE fixed_profit;

scalar sma;
sma=smin(iter,profit(iter));
File TestingFile3 / Alg.csv /;
TestingFile3.pc=5;
TestingFile3.nd=5;
put TestingFile3; 
put 'Omega', put 'Tolerance', put 'Step Size Rule', put 'Iterations', put 'Converged?', put 'Gap LR', put 'Gap Naive', put 'Obj. Naive', put 'Obj. LR', put 'Gap' put 'Time Naive', put 'Time LR', put 'Final Lambda', put 'LB Heuristic' put /;
put '', put '', put '', put '', put '', put '', put '', put '', put '', put '', put '', put '', put '', put sma put /;


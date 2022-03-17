$ONTEXT
This is Step 4 of Algorithm 1 of file v10.pdf (iEEE paper)
Big M is updated
Single file to run different values of rho
Promise is updated at each iteration
Can update with new wind scenarios file

After this file, we run fixed scenario problem - dynamic_1500.gms

$OFFTEXT

OPTIONS PROFILE =3, RESLIM   = 2100, LIMROW   = 5, LP = CPLEX, MIP = gurobi, RMIP=Gurobi, NLP = CONOPT, MINLP = DICOPT, MIQCP = CPLEX,
         SOLPRINT = OFF, decimals = 8, optcr=0.0005, optca=0.0005, threads =8, integer4=1;

********************************************************************************
* run_time_total recorder
SCALAR start_time, end_time, run_time_total;

*-------------------------------------------------------------------------------

SETS T times/t1*t24/;
SETS W scenarios /scen1*scen4272/;
SETS iter iterations /iter1*iter30/;
*SETS rho_ind /r1*r5/;
*parameter rho_val(rho_ind) / r1 5, r2 40, r3 80, r4 120, r5 160/;
*SETS rho_ind /r1*r3/;
*parameter rho_val(rho_ind) / r1 40, r2 50, r3 60/;

SETS rho_ind /r1/;
parameter rho_val(rho_ind) / r1 40/;

ALIAS (T,TT);
ALIAS (W,I);
ALIAS (W,Scen);
** define generator costs and wind selling prices


TABLE Solar(scen,t)
$ondelim
*$INCLUDE %SOLAR%.csv
$INCLUDE ScenarioSet_hierarchical_solar7200.csv
$offdelim
;

*Tolerance 
scalar tol;
*tol=%TOL%;
tol=0.07;


* time limit for each problem
scalar time_limit;
*time_limit=%TIMELIM%;
time_limit=2250;

ALIAS (T,TT);
alias(scen,i);

scalar n;
n=card(scen);

*Scalar which tells if LR converges
scalar convergence;

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


table probability(w,*)
$ONDELIM
$INCLUDE ScenarioWeights_hierarchical_solar7200.csv
$OFFDELIM
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
scalar it ;
set dummy_set(scen);
* make the dum_iter go till at least the size of it
set dum_iter /dum_iter1*dum_iter100/;
scalar G upper bound on q - p ;
scalar run_time_total, start_time, end_time;

scalars m,n,size, step, delta, rho, new_n, new_m, tot_time ;    ;
parameters time(rho_ind,iter), promise(rho_ind,iter,t), profit(rho_ind,iter) , num_scen(rho_ind,iter), y_previous(t) ;

* number of initial scenarios and % increase at each iterations
step = 20;
delta = 0.1 ;



********************************************************************************
*                                begin model
********************************************************************************



POSITIVE VARIABLES P(scen,t), Q(scen,t), Y(T), X(scen,t), dummy(t) ;
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
        
        Const_dum1(t)
        Const_dum2(t)
        ;

Objective.. OBJ=E= SUM(T,Prices(T, 'REW')*Y(T) - Sum(w$((ord(w) ge n) and (ord(w) le m)),probability(w,'value')* ( Prices(T, 'CHAR')* P(w,t) + Prices(t, 'DISCHAR') * Q(w,t) ) ) )  + rho*sum(t,dummy(t))   ;

Const1(scen,t)$(ord(t) lt card(t) and (ord(scen) ge n) and (ord(scen) le m))..
         X(scen,t+1) =E= X(scen,t) + eta* P(scen,t) - (1/eta)* Q(scen,t) ;

Const_chance_1(scen,t)$((ord(scen) ge n) and (ord(scen) le m)).. Y(T) + P(scen,t) -  Q(scen,t) -SOLAR(scen,t) =L= Z(scen)*BigM(scen,t) ;


Const_chance_2..      - sum(scen$((ord(scen) ge n) and (ord(scen) le m)), z(scen)*probability(scen,'value')) =G= -tol;

Const_dum1(t).. dummy(t)=g= y(t) - y_previous(t);
Const_dum2(t).. dummy(t)=g= y_previous(t) - y(t);

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


********************************************************************************
*                                begin iterations
********************************************************************************

*rho_ind -> r1 40
loop(rho_ind,   
tot_time = 0;
size =0;
*We take the first 20 scenarios 
new_n = 1;
new_m = 20;
y_previous(t) = 0;

*iter --> iter1 * iter1000
loop(iter,

if (tot_time < 1800,

*rho_val(rho_ind=r1)=40
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

maxsolar(t) =smax(scen,solar(scen,t)) ;
dummysolar(scen,t) = solar(scen,t) ;

it = floor(card(scen)*tol) + 1;

* index of it
loop(t,
loop(dum_iter$(ord(dum_iter)le it),
* find the smallest solar value for this t
         minsolar(t) = smin(w$((ord(w) ge n) and (ord(w) le m)),dummysolar(w,t)) ;
* index of smallest solar value
         dummy_set(w)$((ord(w) ge n) and (ord(w) le m)) = yes$(dummysolar(w,t) eq minsolar(t)) ;
* make the smallest solar value large
         dummysolar(w,t)$(dummy_set(w) and ((ord(w) ge n) and (ord(w) le m))) =maxsolar(t) ;
); );
G = min(eta*(BigX - LowX), max_discharge) ;

BigM(scen,t)$((ord(scen) ge n) and (ord(scen) le m))= G - solar(scen,t) + minsolar(t);


         start_time=jnow;
         SOLVE SCHEDULE USING MIP MINIMIZING OBJ ;
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

FILE sampled_dynamic /sampled_dynamic_hi.csv/;
sampled_dynamic.PC = 5;
sampled_dynamic.ND = 3;
PUT sampled_dynamic;
loop(rho_ind,LOOP(t, loop(iter, put promise(rho_ind,iter,t);) put /; ); );

PUTCLOSE sampled_dynamic;


FILE profit_value /profit_value_hi.csv/;
profit_value.PC = 5;
profit_value.ND = 3;
PUT profit_value;
loop(rho_ind, put rho_val(rho_ind) put /;
put ''  put 'time' put ' number scenarios' put 'profit_value' put /;
LOOP(iter, put iter.tl PUT time(rho_ind,iter) put num_scen(rho_ind,iter) put profit(rho_ind,iter) put /;);
);
PUTCLOSE profit_value;



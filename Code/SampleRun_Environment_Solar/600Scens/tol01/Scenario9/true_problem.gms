$ONTEXT
March 2019
True problem

NEED TO MANUALLY SORT THE EXCEL FILE!!!!

$OFFTEXT

OPTIONS PROFILE =3, RESLIM   = 2100, LIMROW   = 5, LP = CPLEX, MIP = cplex, RMIP=cplex, NLP = CONOPT, MINLP = DICOPT, MIQCP = CPLEX, SOLPRINT = OFF, decimals = 8, optcr=0.001, optca=0.001, threads =8, integer4=0;

********************************************************************************

** sets later to be defined in input file
SETS T times/t1*t24/;
SETS SCEN scenarios /scen1*scen150/;

ALIAS (T,TT);

** define battery  operation costs costs and solar selling prices

TABLE PRICES(t,*)
$ONDELIM
$INCLUDE battery_revenue.csv
$OFFDELIM
;

** define solar scenarios at all time periods
TABLE Solar(scen,t)
$ondelim
$INCLUDE solar_scenarios_150.csv
$offdelim


* Scaling of Solar power scenarios ;
scalar scale ;
scale = 1;
Solar(scen,t) = scale* Solar(scen,t) ;
* Remove too many decimals in Solar
Solar(scen,t) = round(Solar(scen,t),2) ;


scalar PROBABILITY;
PROBABILITY = 1/CARD(scen);
;
scalar eta ;
*from Ben paper
eta = 0.9
;

parameters max_store(t), min_store(t), max_charge, max_discharge;


** define tolerance threshold
SCALAR tol, threshold;
TOL       = 0.05;
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
set dummy(scen);
* make the dum_iter go till at least the size of it
set dum_iter /dum_iter1*dum_iter100/;
loop(t,
loop(dum_iter$(ord(dum_iter)le it),
* find the smallest solar value for this t
         minsolar(t) = smin(scen,dummysolar(scen,t)) ;
* index of smallest solar value
         dummy(scen) = yes$(dummysolar(scen,t) eq minsolar(t)) ;
* make the smallest solar value large
         dummysolar(scen,t)$dummy(scen) =maxsolar(t) ;
); );
scalar G upper bound on q - p ;
G = min(BigX - LowX, max_discharge)  ;

BigM(scen,t)= G - solar(scen,t) + minsolar(t);
**********************************

********************************************************************************
* Lagrangian dual
* Let Const_chance_2 be the complicating constraint
*---------------------------------------------------------------------
********************************************************************************

parameter lambda;
parameter ldual value of Lagrangian dual ;
parameter bound total value of Lagrangian dual ;
variable bound_lr objective of Lagrangian dual;
scalar init_lambda, init_bound initial value of lambda dual LP objective from LP ;


********************************************************************************
*                                begin model
********************************************************************************
POSITIVE VARIABLES P(scen,t), Q(scen,t), Y(T), X(scen,t) ;
VARIABLES OBJ, BOUND_LR;
BINARY VARIABLE W(scen,t), Z(scen) ;

EQUATIONS
        Objective
        Const1_1(scen,t)    balance constraint
        Const1_2(scen,t)    balance constraint
        Const1(scen,t)
        Const2(scen,t)    max charge
        Const3(scen,t)    max discharge

        ;

Objective.. OBJ=E= SUM(T,Prices(T, 'REW')*Y(T) - probability*Sum(scen, Prices(T, 'CHAR')* P(scen,t) + Prices(t, 'DISCHAR') * Q(scen,t)  ) )    ;

Const1_1(scen,t)$(ord(t) lt card(t)).. LowX =L= X_0 + eta* sum(tt$(ord(tt) le ord(t)),  P(scen,tt)) - (1/eta)*sum(tt$(ord(tt) le ord(t)), Q(scen,tt)) ;
Const1_2(scen,t)$(ord(t) lt card(t)).. BigX =G= X_0 + eta* sum(tt$(ord(tt) le ord(t)),  P(scen,tt)) - (1/eta)*sum(tt$(ord(tt) le ord(t)), Q(scen,tt)) ;
Const1(scen,t)$(ord(t) le card(t))..   X(scen,t+1) =E= X(scen,t) + eta* P(scen,t) - (1/eta)* Q(scen,t) ;

Const2(scen,t).. P(scen,t) =L= max_charge * W(scen,t)  ;
Const3(scen,t).. Q(scen,t) =L= max_discharge * ( 1- W(scen,t) )  ;

*** bounds on any variables
x.up(scen,t) = BigX ;
x.lo(scen,t) = LowX ;
q.up(scen,t) = max_discharge ;
p.up(scen,t) =  max_charge ;
x.fx(scen,'t1') = X_0 ;
z.prior(scen)   = 1;



model schedule     /all/ ;
solve schedule using MIP maximizing Obj ;

display Obj.l;


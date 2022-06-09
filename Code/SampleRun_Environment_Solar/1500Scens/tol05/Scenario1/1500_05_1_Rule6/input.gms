** INPUT SPECIFICATIONS **

SETS T times/t1*t24/;
SETS SCEN scenarios /scen1*scen1500/;

TABLE Solar(scen,t)
$ondelim
$INCLUDE solar_scenarios.csv
$offdelim
;

** define tolerance  
scalar tol;
tol=0.05;

** set maximum number of iterations in LR
set iter number of subgradient iterations /iter1*iter10/;

* Import the SORTED file for the QP bound
table scenario_sorted(scen,*)
$ondelim
$INCLUDE scenario_sorted.csv
$offdelim
;

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

* Scaling of Solar power scenarios ;
scalar scale ;
scale = 1;
Solar(scen,t) = scale* Solar(scen,t) ;
* Remove too many decimals in Solar
Solar(scen,t) = round(Solar(scen,t),2) ;


scalar PROBABILITY;
PROBABILITY = 1/CARD(scen);

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
G = min(eta*(BigX - LowX), max_discharge) ;

BigM(scen,t)= G - solar(scen,t) + minsolar(t);

scalar run_time_total, start_time, end_time, LP_time, bound_time, lr_time ;

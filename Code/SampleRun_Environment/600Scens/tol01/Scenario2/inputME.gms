** sets later to be defined in input file

** to be changed

SETS T times/t1*t24/;
* Number of scenarios 
*SETS SCEN scenarios /scen1*%MAXSCEN%/;
SETS SCEN scenarios /scen1*scen600/;
alias(scen, w);
Parameter last_z(scen);

SCALAR G, GG, start_cost, ramp, uptime, downtime;
G=130;
GG=20;
ramp=50;
UPTIME =3;
DOWNTIME=3;

TABLE wind(scen,t)
$ondelim
*$INCLUDE %SOLAR%.csv
$INCLUDE wind_scenarios_600_2.csv
$offdelim
;
*alias(solar,wind);

*Tolerance 
scalar tol;
*tol=%TOL%;
tol=0.05;

* maximum number of iterations in LR
set iter number of subgradient iterations /iter1*iter10/;

* time limit for each problem
scalar time_limit;
*time_limit=%TIMELIM%;
time_limit=2250;


* Import the SORTED file
$ONTEXT
table scenario_sorted(scen,*)
$ondelim
$INCLUDE scenario_sorted.csv 
$offdelim
;
$OFFTEXT




ALIAS (T,TT);
alias(scen,i);

scalar n;
n=card(scen);

*Scalar which tells if LR converges
scalar convergence;

** define battery  operation costs costs and solar selling prices

TABLE PRICES(t,*)
$ONDELIM
$INCLUDE wind_costs.csv
$OFFDELIM
;
*alias(costs, prices);

Prices(t,'rew')     =  - Prices(t,'rew');
Prices(t,'gen')    =  - Prices(t,'gen');
** define solar scenarios at all time periods


* Scaling of wind power scenarios ;
scalar scale ;
scale = 1;
wind(scen,t) = scale* wind(scen,t) ;
* Remove too many decimals in Solar
wind(scen,t) = round(wind(scen,t),2) ;


scalar PROBABILITY;
PROBABILITY = 1/CARD(scen);


parameters max_store(t), min_store(t), max_charge, max_discharge;


** define tolerance threshold
SCALAR threshold;
threshold = floor(card(scen)*TOL)  ;

$ontext
parameter BigX, LowX, X_0 maximum minimum initial energy stored ;
parameter BigM(scen,t) find a good BigM ;

BigX = 960 ;
LowX = 0.2* BigX ;
X_0  = 0.5* BigX ;
max_charge =  0.5* BigX ;
max_discharge =  0.5* BigX ;
$offtext

************** Find a Big M

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



parameter BigMM(w), BigM(w,t);
BigMM(w)= smax(t, wind(w,t));
BigM(w,t)= G - wind(w,t) + minwind(t);


scalar run_time_total, start_time, end_time, LP_time, bound_time, lr_time ;

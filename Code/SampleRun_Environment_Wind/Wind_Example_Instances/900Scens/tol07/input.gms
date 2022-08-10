** INPUT SPECIFICATIONS **

SETS T times/t1*t24/;
SETS W scenarios /scen1*scen900/;

ALIAS (T,TT);
ALIAS (W,I);
ALIAS (W,SCEN);

** define generator costs and wind selling prices
TABLE PRICES(T,*)
$ONDELIM
$INCLUDE wind_costs.csv
$OFFDELIM
;

Prices(t,'rew')     =  - Prices(t,'rew');
Prices(t,'gen')    =  - Prices(t,'gen');

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


Parameter last_z(scen);

SCALAR G, GG, start_cost, ramp, uptime, downtime;
G=130;
GG=20;
ramp=50;
UPTIME =3;
DOWNTIME=3;

* maximum number of iterations in LR
set iter number of subgradient iterations /iter1*iter10/;
* maximum number of iteration in fixed (bound IR)
set iterFIX iterations for fixed /iter1*iter30/;


* time limit for each problem
scalar time_limit;
time_limit=2250;

scalar n;
n=card(scen);

*Scalar which tells if LR converges
scalar convergence;

* Scaling of wind power scenarios ;
scalar scale ;
scale = 1;
wind(scen,t) = scale* wind(scen,t) ;
* Remove too many decimals in Solar
wind(scen,t) = round(wind(scen,t),2) ;


parameters max_store(t), min_store(t), max_charge, max_discharge;


** define tolerance threshold
SCALAR threshold;
threshold = floor(card(scen)*TOL)  ;


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

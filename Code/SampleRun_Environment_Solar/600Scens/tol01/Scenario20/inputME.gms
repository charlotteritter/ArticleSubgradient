** sets later to be defined in input file

** to be changed

SETS T times/t1*t24/;

* Number of scenarios 
*SETS SCEN scenarios /scen1*%MAXSCEN%/;
SETS SCEN scenarios /scen1*scen600/;

TABLE Solar(scen,t)
$ondelim
*$INCLUDE %SOLAR%.csv
$INCLUDE solar_scenarios_600_20.csv
$offdelim
;

*Tolerance 
scalar tol;
*tol=%TOL%;
tol=0.01;

* maximum number of iterations in LR
set iter number of subgradient iterations /iter1*iter10/;

* time limit for each problem
scalar time_limit;
*time_limit=%TIMELIM%;
time_limit=2250;


* Import the SORTED file
table scenario_sorted(scen,*)
$ondelim
*$INCLUDE %SORTEDFILE%.csv
$INCLUDE scenario_sorted_600_20_01.csv
*$INCLUDE scenario_sorted_100_2_01.csv
*$INCLUDE scenario_sorted_100_3_01.csv
*$INCLUDE scenario_sorted_100_4_01.csv
*$INCLUDE scenario_sorted_100_5_01.csv
*$INCLUDE scenario_sorted_100_6_01.csv
*$INCLUDE scenario_sorted_100_7_01.csv
*$INCLUDE scenario_sorted_100_8_01.csv
*$INCLUDE scenario_sorted_100_9_01.csv
*$INCLUDE scenario_sorted_100_10_01.csv
*$INCLUDE scenario_sorted_100_11_01.csv
*$INCLUDE scenario_sorted_100_12_01.csv
*$INCLUDE scenario_sorted_100_13_01.csv
*$INCLUDE scenario_sorted_100_14_01.csv
*$INCLUDE scenario_sorted_100_15_01.csv
*$INCLUDE scenario_sorted_100_16_01.csv
*$INCLUDE scenario_sorted_100_17_01.csv
*$INCLUDE scenario_sorted_100_18_01.csv
*$INCLUDE scenario_sorted_100_19_01.csv
*$INCLUDE scenario_sorted_100_20_01.csv
*$INCLUDE scenario_sorted_100_1_05.csv
*$INCLUDE scenario_sorted_100_2_05.csv
*$INCLUDE scenario_sorted_100_3_05.csv
*$INCLUDE scenario_sorted_100_4_05.csv
*$INCLUDE scenario_sorted_100_5_05.csv
*$INCLUDE scenario_sorted_100_6_05.csv
*$INCLUDE scenario_sorted_100_7_05.csv
*$INCLUDE scenario_sorted_100_8_05.csv
*$INCLUDE scenario_sorted_100_9_05.csv
*$INCLUDE scenario_sorted_100_10_05.csv
*$INCLUDE scenario_sorted_100_11_05.csv
*$INCLUDE scenario_sorted_100_12_05.csv
*$INCLUDE scenario_sorted_100_13_05.csv
*$INCLUDE scenario_sorted_100_14_05.csv
*$INCLUDE scenario_sorted_100_15_05.csv
*$INCLUDE scenario_sorted_100_16_05.csv
*$INCLUDE scenario_sorted_100_17_05.csv
*$INCLUDE scenario_sorted_100_18_05.csv
*$INCLUDE scenario_sorted_100_19_05.csv
*$INCLUDE scenario_sorted_100_20_05.csv
*$INCLUDE scenario_sorted_150_01.csv
*$INCLUDE scenario_sorted_150_03.csv
*$INCLUDE scenario_sorted_150_05.csv
*$INCLUDE scenario_sorted_300_01.csv
*$INCLUDE scenario_sorted_300_03.csv
*$INCLUDE scenario_sorted_300_05.csv
*$INCLUDE scenario_sorted_450_01.csv
*$INCLUDE scenario_sorted_450_03.csv
*$INCLUDE scenario_sorted_450_05.csv
*$INCLUDE scenario_sorted_600_01.csv
*$INCLUDE scenario_sorted_600_05.csv
*$INCLUDE scenario_sorted_600_03.csv
*$INCLUDE scenario_sorted_900_01.csv
*$INCLUDE scenario_sorted_900_03.csv
*$INCLUDE scenario_sorted_900_05.csv
*$INCLUDE scenario_sorted_1200_01.csv
*$INCLUDE scenario_sorted_1200_03.csv
*$INCLUDE scenario_sorted_1200_05.csv
*$INCLUDE scenario_sorted_2400_01.csv
*$INCLUDE scenario_sorted_2400_05.csv
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
** define solar scenarios at all time periods


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

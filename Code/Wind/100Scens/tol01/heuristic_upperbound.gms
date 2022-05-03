********************************************************************************
* Find a lower bound on the problem : a feasible solution
********************************************************************************
Prices(t,'rew')     =  - Prices(t,'rew');
Prices(t,'gen')    =  - Prices(t,'gen');

* x values for the fixed model (IR)
TABLE y_100(t,iterFIX)
$ONDELIM
$INCLUDE sampled_dynamic.csv
$OFFDELIM
;

start_time = jnow ;
schedule_scenario.solvelink = 5 ;
alias(rs,scen)

*********************************************************************************
* Solve the fixed problem of the SAA to get an upperbound to the problem
*********************************************************************************

scalar temp;
scalar currentMAX;
currentMAX=0;

loop(iterFIX,
         y.fx(t)=y_100(t,iterFIX);
         t1=jnow ;
         SOLVE SCHEDULEFIX USING MIP MAXIMIZING OBJ;
         t2=jnow;
         run_time(iterFIX) = ghour(t2 - t1)*3600 + gminute(t2 - t1)*60 + gsecond(t2 - t1);
         profit(iterFIX)= obj.l;
         temp=obj.l;
         if(temp>currentMAX, currentMAX=temp; prev_y(t) = y.l(t); prev_z(scen)=z.l(scen););

display prev_y,prev_z;
);

********************************************************************************
*                                write output
********************************************************************************
ALIAS(iterFIX, iterF);
profitFIXED=smax(iterF, profit(iterF));
upperbound=-profitFIXED;

Prices(t,'rew')     =  - Prices(t,'rew');
Prices(t,'gen')    =  - Prices(t,'gen');



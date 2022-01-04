********************************************************************************
* Find a upperbound on the problem : a feasible solution
********************************************************************************
Prices(t,'rew')     =  - Prices(t,'rew');
Prices(t,'gen')    =  - Prices(t,'gen');


*upperbound =  0;
* Find a upper bound using a fixed value and solving MIP (a feasible solution)
* solve single scenario problem and choose the worst #threshold problems
start_time = jnow ;
schedule_scenario.solvelink = 5 ;
alias(rs,scen)

*********************************************************************************
* Solve the fixed problem of the hierarchical SAA to get an upperbound to the problem
*********************************************************************************


* 1st solve with the plain old 20 scenario promise
* use a solution from 20 scenarios scenarios
*y.fx(t)=y_20(t,'values');
*SOLVE SCHEDULE USING MIP MAXIMIZING OBJ;
*profit_orig= obj.l;
scalar temp;
scalar currentMAX;
currentMAX=0;
*2nd iteration onward

loop(iterFIX,
         y.fx(t)=y_100(t,iterFIX);
         t1=jnow ;
         SOLVE SCHEDULEFIX USING MIP MAXIMIZING OBJ;
         t2=jnow;
         run_time(iterFIX) = ghour(t2 - t1)*3600 + gminute(t2 - t1)*60 + gsecond(t2 - t1);
         profit(iterFIX)= obj.l;
*        profit(iterFIX)=10;
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



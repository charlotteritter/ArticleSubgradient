********************************************************************************
* Find a lower bound on the problem : a feasible solution (QP)
********************************************************************************

* Algorithm 4
* Find a upper bound using a fixed value and solving MIP (a feasible solution)
* solve single scenario problem and choose the worst #threshold problems
start_time = jnow ;
schedule_scenario.solvelink = 5 ;
alias(rs,scen)
parameter res_scenarios(rs,*) ;

loop(rs,
counter = ord(rs)  ;
*Solving (2) for each scenario w with corresponding z(w) = 0
solve schedule_scenario using MIP minimizing Obj ;
res_scenarios(rs,'obj') =obj.l ;


res_scenarios(rs,'scenario') = counter ;

);
end_time =jnow ;
run_time_total = ghour(end_time - start_time)*3600 + gminute(end_time - start_time)*60 + gsecond(end_time - start_time);


** Write to a file to sort
FILE scen_sorted /scen_sorted.csv/;
*.PC specifies the format of the put file: 5 -> Formatted output; non-numeric output is quoted and each item is delimited with commas.
scen_sorted.PC = 5;
*.ND sets the number of decimals that are displayed in the put file
scen_sorted.ND = 3;
PUT scen_sorted;
loop(rs, put res_scenarios(rs,'scenario') put res_scenarios(rs,'obj') put /; ) ;
PUTCLOSE scen_sorted;


z.fx(scen) = scenario_sorted(scen,'value') ;
*** Ensure file was generated correctly
if ( sum(scen,z.lo(scen)) ne threshold, abort "sorted file not generated correctly check manually") ;
start_time = jnow ;
schedule.solprint = 0;
*schedule.optfile  = 1;
schedule.solvelink = 5 ;
solve schedule using MIP minimizing Obj ;
end_time = jnow ;
bound_time =  run_time_total + ghour(end_time - start_time)*3600 + gminute(end_time - start_time)*60 + gsecond(end_time - start_time);
upperbound = Obj.l ;
prev_y(t) = y.l(t) ;
* Clear bound on z now
z.up(scen) = 1 ;
z.lo(scen) = 0 ;

display upperbound,prev_y  ;

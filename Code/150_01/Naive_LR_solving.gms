$ONTEXT
Bismark Singh
March 17, 2021

Code for plain Lagrangian decomposition (removed Progressive Hedging)
Based on the paper at http://www.optimization-online.org/DB_HTML/2019/05/7222.html
 

Excel file used for LB heuristic needs to be manually sorted

$OFFTEXT

$eolcom //
OPTIONS PROFILE =3, RESLIM   = 4200, LIMROW   = 5, LP = CPLEX, MIP = gurobi, RMIP=Gurobi, NLP = CONOPT, MINLP = DICOPT, MIQCP = CPLEX,SOLPRINT = OFF, decimals = 8, optcr=0.01, optca=0.01, threads =8, integer4=0;

********************************************************************************
*                                Include input files
********************************************************************************
$include inputME.gms // no need to change for Lagrangian decomposition
$include subgradient_parameters.gms

$include equations_all.gms


scalar d;
*set indices /1/;
set indices /1*6/;

File TestingFile3 / TestingFile3.csv /;
TestingFile3.pc=5;
TestingFile3.nd=5;
put TestingFile3; 
put 'Omega', put 'Tolerance', put 'Step Size Rule', put 'Iterations', put 'Converged?', put 'Gap LR', put 'Gap Naive', put 'Obj. Naive', put 'Obj. LR', put 'Gap' put 'Time Naive', put 'Time LR', put 'Final Lambda', put 'LB Heuristic' put /;

********************************************************************************
* Solve main Problem
********************************************************************************
option limrow = 10000;
schedule.optfile=1;

start_time = jnow;
solve schedule using MIP minimizing Obj ;                       
end_time = jnow ;



run_time_total = ghour(end_time - start_time)*3600 + gminute(end_time - start_time)*60 + gsecond(end_time - start_time);

scalar ObjNaive;
ObjNaive=Obj.l;

scalar zlower;
zlower=-Obj.l;

scalar zupper;
zupper=-schedule.objEst;

scalar GapNaive;
GapNaive = (zupper-zlower)/zupper;

scalar ObjLR;

scalar heuristic;

scalar TimeNaive;
TimeNaive=run_time_total;    

display Obj.l, run_time_total ;

********************************************************************************
* Solve the Lagrangian Dual problem now
********************************************************************************

$include lp_lowerbound.gms // no need to change for Lagrangian decomposition
$include heuristic_upperbound.gms // no need to change for Lagrangian decomposition

parameter ldual_iter(iter) obj function at each iteration ;
lr_time = 0 ;

option limrow = 0, limcol = 0, optca=0.0001, optcr=0.0001, RESLIM   = 2100;

prev_y(t) = y.l(t) ;

parameter check(scen,t);
scalar steprule;
scalar FinalIter;

loop(indices,
    option clear=results;
    noimprovement = 0;
    lambda=init_lambda;
    lowerbound=LP_bound;
    theta=originalTheta;
    lr_time=0;
    run_time_total=0;
    contin=1;
    steprule=ord(indices);
    
    loop(iter$contin,
    num_iter = ord(iter) ;
*         pass a warm start
             y.l(t) = prev_y(t) ;
             z.l(scen) = scenario_sorted(scen,'value') ;
             start_time = jnow;
    
*********************************************************************
***Solve a Lagrangian iteration 
*********************************************************************

*Test
    
$include plain_lr.gms
    
    end_time = jnow ;
    results(iter,'time') = ghour(end_time - start_time)*3600 + gminute(end_time - start_time)*60 + gsecond(end_time - start_time);
    results(iter,'objective') = bound ;
    
$include LR_updates.gms
    if( ((results(iter,'gap') < exit_tol) and (num_iter > 2)),convergence=2; contin = 0;);
    lr_time = lr_time + results(iter,'time')   ;
    if (lr_time > 2400, contin = 0 ; ) ;
    
    d=results(iter,'gap');
    FinalIter=num_iter;
);
    
run_time_total = LP_time + lr_time + bound_time  ;
    
* check if any p and q active simultaneously (nothing to do with Lagrangian)
*parameter check(scen,t);


ObjLR=-lowerbound;
heuristic=-upperbound;

put TestingFile3;
put n, put tol, put steprule, put FinalIter, put convergence, put d, put GapNaive, put zlower, put ObjLR, put ((ObjLR-max(heuristic,zlower))/ObjLR), put TimeNaive, put lr_time, put lambda, put heuristic put /;

display results, lowerbound, upperbound, LP_bound, run_time_total, lr_time, num_iter ;
display z.l, y.l ;
display zlower, zupper, ObjLR, heuristic;

);


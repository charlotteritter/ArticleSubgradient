ONTEXT
Bismark Singh
March 17, 2021

Code for plain Lagrangian decomposition (removed Progressive Hedging)
Based on the paper at http://www.optimization-online.org/DB_HTML/2019/05/7222.html
 

Excel file used for LB heuristic needs to be manually sorted

$OFFTEXT

$eolcom //
OPTIONS PROFILE =3, RESLIM   = 2100, LIMROW   = 5, LP = CPLEX, MIP = cplex, RMIP=cplex, NLP = CONOPT, MINLP = DICOPT, MIQCP = CPLEX, SOLPRINT = OFF, decimals = 8, optcr=0.001, optca=0.001, threads =8, integer4=0;

********************************************************************************
*                                Include input files
********************************************************************************
$include inputME.gms // no need to change for Lagrangian decomposition
$include subgradient_parameters.gms

$include equations_all.gms
$include lp_lowerbound.gms // no need to change for Lagrangian decomposition
$include heuristic_upperbound.gms // no need to change for Lagrangian decomposition

scalar r;
set indices /1*6/;

File TestingFile / TestingFile.csv /;
TestingFile.pc=5;
TestingFile.nd=5;
put TestingFile; 
put 'Omega', put 'Tolerance', put 'Step Size Rule', put 'Iterations', put 'Converged?', put 'Gap LR', put 'Gap Naive', put 'Obj. Naive', put 'Obj. LR', put 'Gap' put 'Time Naive', put 'Time LR' put 'Final Lambda' put/;

********************************************************************************
* Solve main Problem
********************************************************************************

start_time = jnow;
solve schedule using MIP minimizing Obj ;
end_time = jnow ;



run_time_total = ghour(end_time - start_time)*3600 + gminute(end_time - start_time)*60 + gsecond(end_time - start_time);

scalar ObjNaive;
ObjNaive=Obj.l;


scalar GapNaive;
GapNaive = abs((schedule.objEst-schedule.objVal)/schedule.objEst);

scalar TimeNaive;
TimeNaive=run_time_total;    

display Obj.l, run_time_total ;

********************************************************************************
* Solve the Lagrangian Dual problem now
********************************************************************************

parameter ldual_iter(iter) obj function at each iteration ;
lr_time = 0 ;

option limrow = 0, limcol = 0, optca=0.0001, optcr=0.0001 ;

prev_y(t) = y.l(t) ;

parameter check(scen,t);
scalar steprule;
scalar FinalIter;

loop(indices,
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
    
$include LR_updatesMe.gms
    if( ((results(iter,'gap') < exit_tol) and (num_iter > 2)),convergence=2; contin = 0;);
    lr_time = lr_time + results(iter,'time')   ;
    if (lr_time > time_limit, contin = 0 ; ) ;
    
    r=results(iter,'gap');
    FinalIter=num_iter;
);
    
run_time_total = LP_time + lr_time + bound_time  ;
    
* check if any p and q active simultaneously (nothing to do with Lagrangian)
*parameter check(scen,t);
check(scen,t) = 0 ;
check(scen,t) = 1$( p.l(scen,t) gt 0 and q.l(scen,t) gt 0) ;
if ( sum((scen,t), check(scen,t)) gt 0, abort "error: p and q are one together, check. ");



put TestingFile;
put n, put tol, put steprule, put FinalIter, put convergence, put r, put GapNaive, put ObjNaive, put lowerbound, put ((lowerbound-ObjNaive)/ObjNaive), put TimeNaive, put lr_time put lambda put /;

display results, lowerbound, upperbound, LP_bound, run_time_total, lr_time, num_iter ;
display z.l, y.l ;


);


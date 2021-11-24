$ONTEXT
Bismark Singh
March 17, 2021

Code for plain Lagrangian decomposition (removed Progressive Hedging)
Based on the paper at http://www.optimization-online.org/DB_HTML/2019/05/7222.html
 

Excel file used for LB heuristic needs to be manually sorted

$OFFTEXT

$eolcom //
OPTIONS PROFILE =3, RESLIM   = 2100, LIMROW   = 5, LP = CPLEX, MIP = gurobi, RMIP=Gurobi, NLP = CONOPT, MINLP = DICOPT, MIQCP = CPLEX, SOLPRINT = OFF, decimals = 8, optcr=0.001, optca=0.001, threads =8, integer4=0;

********************************************************************************
*                                Include input files
********************************************************************************
$ include inputME.gms // no need to change for Lagrangian decomposition
$include subgradient_parameters.gms

$include equations_all.gms
$include lp_lowerbound.gms // no need to change for Lagrangian decomposition
$include heuristic_upperbound.gms // no need to change for Lagrangian decomposition

File TestingFile3 / LR5.csv /;
TestingFile3.pc=5;
TestingFile3.nd=5;
put TestingFile3; 
put 'Omega', put 'Tolerance', put 'Step Size Rule', put 'Iterations', put 'Converged?', put 'Gap LR', put 'Gap Naive', put 'Obj. Naive', put 'Obj. LR', put 'Gap' put 'Time Naive', put 'Time LR', put 'Final Lambda', put 'LB Heuristic' put /;

********************************************************************************
* Solve the Lagrangian Dual problem now
********************************************************************************
scalar d;
scalar FinalIter;
sets ro "rows" /1/;
sets cl "columns" /1*3/;
$call csv2gdx Naive.csv id=nv index=1 values=2..lastCol useHeader=y
$gdxIn Naive.gdx
*$load i=dim1
*$load j=dim2
parameter nv(ro,cl) ;
$load nv
$gdxIn
display ro,cl,nv;

scalar zlower;
zlower=nv('1','1');

scalar GapNaive;
GapNaive =nv('1','2');

scalar TimeNaive;
TimeNaive= nv('1','3');

parameter ldual_iter(iter) obj function at each iteration ;
lr_time = 0 ;

option limrow = 0, limcol = 0, optca=0.0001, optcr=0.0001 ;

prev_y(t) = y.l(t) ;
scalar steprule;
steprule=5;
loop(iter$contin,
num_iter = ord(iter) ;
*         pass a warm start
         y.l(t) = prev_y(t) ;
         z.l(scen) = scenario_sorted(scen,'value') ;
         start_time = jnow;

*********************************************************************
***Solve a Lagrangian iteration 
*********************************************************************
$include plain_lr.gms

         end_time = jnow ;
         results(iter,'time') = ghour(end_time - start_time)*3600 + gminute(end_time - start_time)*60 + gsecond(end_time - start_time);
         results(iter,'objective') = bound ;

$include LR_updates.gms
         if( ((results(iter,'gap') < 0.001) and (num_iter > 2)), contin = 0;);
         lr_time = lr_time + results(iter,'time')   ;
         if (lr_time > 2250, contin = 0 ;) ;
         d=results(iter,'gap');
         FinalIter=num_iter;

);

run_time_total = LP_time + lr_time + bound_time  ;
scalar ObjLR;
scalar heuristic;

ObjLR=-lowerbound;
heuristic=-upperbound;

display results, lowerbound, upperbound, LP_bound, run_time_total, lr_time, num_iter ;
display z.l, y.l ;
display ObjLR, heuristic;

put TestingFile3;
put n, put tol, put steprule, put FinalIter, put convergence, put d, put GapNaive, put zlower, put ObjLR, put ((ObjLR-max(heuristic,zlower))/ObjLR), put TimeNaive, put lr_time, put lambda, put heuristic put /;

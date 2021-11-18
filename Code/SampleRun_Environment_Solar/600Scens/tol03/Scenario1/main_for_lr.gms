$ONTEXT
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
$include input.gms // no need to change for Lagrangian decomposition
$include subgradient_parameters.gms

$include equations_all.gms
$include lp_lowerbound.gms // no need to change for Lagrangian decomposition
$include heuristic_upperbound.gms // no need to change for Lagrangian decomposition

********************************************************************************
* Solve the Lagrangian Dual problem now
********************************************************************************

scalar steprule /6/;

parameter ldual_iter(iter) obj function at each iteration ;
lr_time = 0 ;

option limrow = 0, limcol = 0, optca=0.0001, optcr=0.0001 ;

prev_y(t) = y.l(t) ;

loop(iter$contin,
num_iter = ord(iter) ;
*         pass a warm start
         y.l(t) = prev_y(t) ;
         z.l(scen) = scenario_sorted(scen,'value') ;
         start_time = jnow;

*********************************************************************
***Solve a Lagrangian iteration 
*********************************************************************

$include plain_LR.gms

         end_time = jnow ;
         results(iter,'time') = ghour(end_time - start_time)*3600 + gminute(end_time - start_time)*60 + gsecond(end_time - start_time);
         results(iter,'objective') = bound ;

*$include LR_updates.gms
$include LR_updatesMe.gms
         if( ((results(iter,'gap') < 0.001) and (num_iter > 2)), contin = 0;);
         lr_time = lr_time + results(iter,'time')   ;
         if (lr_time > 2250, contin = 0 ; ) ;
);

run_time_total = LP_time + lr_time + bound_time  ;

* check if any p and q active simultaneously (nothing to do with Lagrangian)
parameter check(scen,t);
check(scen,t) = 0 ;
check(scen,t) = 1$( p.l(scen,t) gt 0 and q.l(scen,t) gt 0) ;
if ( sum((scen,t), check(scen,t)) gt 0, abort "error: p and q are one together, check. ")

display results, lowerbound, upperbound, LP_bound, run_time_total, lr_time, num_iter ;
display z.l, y.l ;


$gdxout results
$unload n tol upperbound steprule lr_time results
execute 'gdxxrw.exe results.gdx var=n'


*execute_unload "results.gdx" n,tol, upperbound, steprule, lr_time, results
*execute 'gdxxrw.exe results.gdx o=results.xls var=x.l'
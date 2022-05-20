$ONTEXT
This is Algorithm 1 with step size 6
$OFFTEXT

OPTIONS PROFILE =3, RESLIM   = 2100, LIMROW   = 5, LP = CPLEX, MIP = cplex, RMIP=cplex, NLP = CONOPT, MINLP = DICOPT, MIQCP = CPLEX, SOLPRINT = OFF, decimals = 8, optcr=0.01, optca=0.01, threads =8, integer4=0;

********************************************************************************
*-------------------------------------------------------------------------------

$include input.gms
$include subgradient_parameters.gms
$include equations_all.gms
$include lp_lowerbound.gms
$include heuristic_upperbound.gms 

*Create output csv file named LR6.csv
File TestingFile3 / LR6.csv /;
TestingFile3.pc=5;
TestingFile3.nd=5;
put TestingFile3; 
put 'Omega', put 'Tolerance', put 'Step Size Rule', put 'Iterations', put 'Converged?', put 'Gap LR', put 'Gap Naive', put 'Obj. Naive', put 'Obj. LR (UB)', put 'Obj. LR (LB),' put 'Gap' put 'Time Naive', put 'Time LR', put 'Final Lambda', put 'LB Heuristic', put 'BBP', put 'Best iter', put 'Lambda best iter', put 'Gamma best iter' /;


********************************************************************************
* Include the Solving of the Naive model 
********************************************************************************
scalar d;
scalar FinalIter;
sets ro "rows" /1/;
sets cl "columns" /1*3/;
$call csv2gdx Naive.csv id=nv index=1 values=2..lastCol useHeader=y
$gdxIn Naive.gdx
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

option limrow = 0, limcol = 0, optca=0.0001, optcr=0.0001;


parameter check(scen,t);
scalar steprule;
steprule=6;
scalar FinalIter;
scalar it;
scalar BBP;
BBP=0;
scalar boundTemp;
scalar lbLR;
scalar gammaBest;
scalar lambdaBest;


    lr_time=0;
    run_time_total=0;
    y.lo(t)=0;
    y.up(t)=INF;
    option clear=y;
    Lagrangian.solveOpt=2;
    loop(iter$contin,
    num_iter = ord(iter) ;
             start_time = jnow;
    
*********************************************************************
***Solve a Lagrangian iteration 
*********************************************************************
    
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

scalar ObjLR;
scalar heuristic;
ObjLR=-lowerbound;
heuristic=-upperbound;
lbLR=-lbLR;

display results, lowerbound, upperbound, LP_bound, run_time_total, lr_time, num_iter ;
display z.l, y.l ;
display zlower, ObjLR, heuristic;


put TestingFile3;
put n, put tol, put steprule, put FinalIter, put convergence, put d, put GapNaive, put zlower, put ObjLR, put lbLR , put ((ObjLR-max(heuristic,zlower))/ObjLR), put TimeNaive, put lr_time, put lambda, put heuristic, put BBP, put it, put lambdaBest, put gammaBest put /;

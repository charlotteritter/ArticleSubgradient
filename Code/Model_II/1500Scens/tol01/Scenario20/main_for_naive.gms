$ONTEXT
This is the GAMS code corresponding to the paper "Statistical performance of subgradient step-size update rules in Lagrangian relaxations of chance-constrained optimization models" by C. Ritter and B. Singh
See: https://github.com/charlotteritter/ArticleSubgradient

This is the naive solution method with a time limit of 4200s for wind model (Model II).
$OFFTEXT

$eolcom //
OPTIONS PROFILE =3, RESLIM   = 4200, LIMROW   = 5, LP = CPLEX, MIP = gurobi, RMIP=gurobi, NLP = CONOPT, MINLP = DICOPT, MIQCP = CPLEX, SOLPRINT = OFF, decimals = 8, optcr=0.001, optca=0.001, threads =8, integer4=0;

********************************************************************************
$include input.gms
$include subgradient_parameters.gms
$include equations.gms   

********************************************************************************
*                                begin model
********************************************************************************


start_time = jnow;
solve schedule using MIP minimizing Obj ;
end_time = jnow ;

run_time_total = ghour(end_time - start_time)*3600 + gminute(end_time - start_time)*60 + gsecond(end_time - start_time);
display Obj.l, run_time_total ;

scalar ObjNaive;
ObjNaive = -Obj.l;

scalar zlower;
zlower=-Obj.l;

scalar zupper;
zupper=-schedule.objEst;

scalar GapNaive;
GapNaive = (zupper-zlower)/zupper;

display ObjNaive;

*Create Output file named Naive.csv, which we include in the LR solving
File TestingFile3 / Naive.csv /;
TestingFile3.pc=5;
TestingFile3.nd=5;
put TestingFile3; 
put '', put "1", put "2", put "3" put /;
put "1", put zlower, put GapNaive, put run_time_total put /;


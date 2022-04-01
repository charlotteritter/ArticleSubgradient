$ONTEXT
Splitting into multiple files

True problem
ELiminnated X variables

$OFFTEXT

$eolcom //
OPTIONS PROFILE =3, RESLIM   = 4200, LIMROW   = 5, LP = CPLEX, MIP = gurobi, RMIP=gurobi, NLP = CONOPT, MINLP = DICOPT, MIQCP = CPLEX, SOLPRINT = OFF, decimals = 8, optcr=0.01, optca=0.01, threads =8, integer4=0;

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


File TestingFile3 / Naive.csv /;
TestingFile3.pc=5;
TestingFile3.nd=5;
put TestingFile3; 
put '', put "1", put "2", put "3" put /;
put "1", put zlower, put GapNaive, put run_time_total put /;

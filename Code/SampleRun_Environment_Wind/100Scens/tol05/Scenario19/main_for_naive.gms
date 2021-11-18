$ONTEXT
Splitting into multiple files

True problem
ELiminnated X variables

$OFFTEXT

$eolcom //
OPTIONS PROFILE =3, RESLIM   = 2100, LIMROW   = 5, LP = CPLEX, MIP = cplex, RMIP=cplex, NLP = CONOPT, MINLP = DICOPT, MIQCP = CPLEX, SOLPRINT = OFF, decimals = 8, optcr=0.001, optca=0.001, threads =8, integer4=0;

********************************************************************************
$ include inputME.gms

** selbst eingefuegt, da sonst error 
$include subgradient_parameters.gms

** vorher: equation_true.gms
$include equations_all.gms   

********************************************************************************
*                                begin model
********************************************************************************
set indices /1*6/;
Set TableSet /1*8/;

TABLE A(indices, TableSet);

start_time = jnow;
solve schedule using MIP minimizing Obj ;
end_time = jnow ;

run_time_total = ghour(end_time - start_time)*3600 + gminute(end_time - start_time)*60 + gsecond(end_time - start_time);

A(indices,'1')=Obj.l;
A(indices,'2')=run_time_total;

execute_unload 'resultsNaive2', A;
execute '=gdx2xls resultsNaive2.gdx';
execute '=shellExecute resultsNaive2.xls';

*execute_unload 'resultsNaive', run_time_total=TimeNaive, Obj.l=ObjNaive;
*execute '=gdx2xls resultsNaive.gdx';
*execute '=shellExecute resultsNaive.xls';

display Obj.l, run_time_total ;




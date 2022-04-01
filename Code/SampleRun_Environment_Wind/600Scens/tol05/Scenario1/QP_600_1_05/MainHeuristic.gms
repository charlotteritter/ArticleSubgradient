$ONTEXT
This is Step 5 of Algorithm 1 of file v10.pdf (iEEE paper)
The 1500 scenario fixed problem
Follow up of SAA.gms

$OFFTEXT

OPTIONS PROFILE =3, RESLIM   = 2100, LIMROW   = 5, LP = CPLEX, MIP = cplex, RMIP=gurobi, NLP = CONOPT, MINLP = DICOPT, MIQCP = CPLEX, SOLPRINT = OFF, decimals = 8, optcr=0.01, optca=0.01, threads =8, integer4=0;

********************************************************************************
*-------------------------------------------------------------------------------

$include input.gms
$include subgradient_parameters.gms
$include equations.gms
*$include lp_lowerbound.gms
$include heuristic_upperbound.gms



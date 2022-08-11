$ONTEXT
This is Step 5 of Algorithm 1 of file v10.pdf (iEEE paper)
The 1500 scenario fixed problem
Follow up of SAA.gms

$OFFTEXT

OPTIONS PROFILE =3, RESLIM   = 1800, LIMROW   = 5, LP = CPLEX, MIP = gurobi, RMIP=Gurobi, NLP = CONOPT, MINLP = DICOPT, MIQCP = CPLEX,
         SOLPRINT = OFF, decimals = 8, optcr=0.0, optca=0.0, threads =8, integer4=1;

********************************************************************************
* Inputs
********************************************************************************

$include input.gms

alias(scen,i);

TABLE y_100(t,iterFIX)
$ONDELIM
$INCLUDE sampled_dynamic.csv
$OFFDELIM
;

* For equations insertion
scalar lambda;

********************************************************************************
*                                begin model
********************************************************************************

$include equations.gms

parameter  profit(iterFIX), y_previous(t), run_time(iterFIX);
scalar profit_orig, t1, t2,start_time,end_time,tot_time;

run_time_total=0;
start_time = jnow;
*2nd iteration onward
loop(iterFIX,
         y.fx(t)=y_100(t,iterFIX);
         t1=jnow ;
         SOLVE SCHEDULE USING MIP MINIMIZING OBJ;
         t2=jnow;
         run_time(iterFIX) = ghour(t2 - t1)*3600 + gminute(t2 - t1)*60 + gsecond(t2 - t1);
         profit(iterFIX)= obj.l;
);
end_time = jnow ;
run_time_total = ghour(end_time - start_time)*3600 + gminute(end_time - start_time)*60 + gsecond(end_time - start_time);
display y.l;

********************************************************************************
*                                write output
********************************************************************************


FILE fixed_profit /fixed_profit.csv/;
fixed_profit.PC = 5;
fixed_profit.ND = 3;
PUT fixed_profit;
loop(iterFIX, put iterFIX.tl put profit(iterFIX) put run_time(iterFIX) put /; );
PUTCLOSE fixed_profit;

scalar sma;
sma=smin(iterFIX,profit(iterFIX));
File TestingFile3 / Alg.csv /;
TestingFile3.pc=5;
TestingFile3.nd=5;
put TestingFile3; 
put 'Omega', put 'Tolerance', put 'Step Size Rule', put 'Iterations', put 'Converged?', put 'Gap LR', put 'Gap Naive', put 'Obj. Naive', put 'Obj. LR', put 'Gap' put 'Time Naive', put 'Time LR', put 'Final Lambda', put 'LB Heuristic' put /;
put '', put '', put '', put '', put '', put '', put '', put '', put '', put '', put '', put '', put run_time_total, put sma put /;



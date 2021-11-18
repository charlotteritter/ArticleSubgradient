scalar tol;
tol=0.05;
scalar omega;
omega=150;
SETS rows /1*150/;
Sets Attributes /scen, value/;


***scenario_sorted was manually created with one column with indices nd descending order of the values
table scenario_sorted(rows,Attributes)
$ondelim
$INCLUDE scenario_sorted.csv
$offdelim
;
scalar c;

scalar Flor;
flor=floor(tol*omega);


***Filling out the zero/one values in the table
loop(rows, c=ord(rows); if(c <= flor, scenario_sorted(rows,'value')=1; else scenario_sorted(rows,'value')=0;););

*display scenario_sorted;

***Creating the .csv
FILE scenario_sorted_X_Y /scenario_sorted_X_Y.csv/;
scenario_sorted_X_Y.PC=5;
scenario_sorted_X_Y.ND=0;
put scenario_sorted_X_Y;
put '', put 'value' put /;
loop(rows, put scenario_sorted(rows, 'scen') put scenario_sorted(rows, 'value') put /;);
putclose scenario_sorted_X_Y;




* Plain LR solve
Lagrangian.solvelink = 5;
Lagrangian.optfile   = 1;
solve Lagrangian using mip minimizing bound_lr  ;

*boundTemp             = bound_lr.l ;
*Save the upper bound of the LR iteration into bound
bound             = Lagrangian.objEst;
display bound ;
prev_y(t)         = y.l(t) ;
last_z(scen)      = z.l(scen) ;

display bound;
* if model is unbounded
if (Lagrangian.modelstat =18,  bound = -100000000;  ); 
results(iter,'status') = Lagrangian.modelstat;
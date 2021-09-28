* Plain LR solve
Lagrangian.solvelink = 5;
Lagrangian.optfile   = 1;
solve Lagrangian using mip minimizing bound_lr  ;

bound             = bound_lr.l ;
prev_y(t)         = y.l(t) ;
last_z(scen)      = z.l(scen) ;

* if model is unbounded
if (Lagrangian.modelstat =18,  bound = -100000000;  ); 
results(iter,'status') = Lagrangian.modelstat;

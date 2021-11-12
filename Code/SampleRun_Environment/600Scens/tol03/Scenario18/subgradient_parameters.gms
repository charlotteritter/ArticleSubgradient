parameter lambda;
parameter ldual value of Lagrangian dual ;
parameter bound total value of Lagrangian dual ;
scalar init_lambda, init_bound initial value of lambda dual LP objective from LP ;
*Hier werden die maximalen Iterationen, also big M festgelegt
*set iter                 number of subgradient iterations /iter1*iter3/; 
scalar num_iter          how many iterations we did ;
scalar contin            stopping             /1/;
parameter stepsize;
scalar theta /2/;
scalar originalTheta;
originalTheta=theta;
scalar noimprovement /0/;
scalar upperbound ;
parameter gamma           subgradient          ;
parameter b;
parameter norm;
parameter norm2;
scalar lowerbound;
parameter lambdaprevious, deltalambda, results(iter,*), prev_w(scen,t), prev_y(t)  ;


scalar m ;
parameter   y_previous(scen,t), y_average_previous(t), weight_previous(scen,t);
scalar profit_orig, t1, t2, exit_tol;

scalar final_gap ;

parameter lb(scen) ;
parameter rho(t) ;

scalar exit_tol ;
exit_tol = 0.00001 ;

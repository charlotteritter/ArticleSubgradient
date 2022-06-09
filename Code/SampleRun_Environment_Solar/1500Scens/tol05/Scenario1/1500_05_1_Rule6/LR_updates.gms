*** UPDATES FOR LAMBDA, GAMMA, THETA

* Gamma update
gamma =  threshold - sum(scen, last_z(scen));
display gamma;

* Bound and theta Update -> in stepsize 6
if (bound > lowerbound,
         lowerbound = bound;
         noimprovement = 0;
         it=num_iter;
         display it;
else
         noimprovement = noimprovement + 1;
         if (noimprovement > 1,  theta = theta/2; noimprovement = 0; );
);

$include stepsizes.gms

* Lambda Update

lambdaprevious = lambda ;

         if (gamma ge 0 and lambdaprevious eq 0,
                 lambda = lambdaprevious ; );
         if (gamma gt 0 and lambdaprevious gt 0,
                 lambda = max(0,lambdaprevious - min(stepsize, lambdaprevious/gamma)*gamma ); );
         if (gamma le 0,
                 lambda = lambdaprevious - stepsize*gamma; );
display lambda;

*Checking if the (modified) BBP holds with equality if the current iteration of Alg. 1 is the "best" one (so the one with the lowest upper bound on the LR until now)
*BBP=0 if BBP doesn't hold with equality in currently best iter, BBP=1 if it does
if (it eq num_iter,
    if (gamma ge 0 and lambdaprevious*gamma eq 0,
        BBP=1;
    else
        BBP=0;
        );
    display BBP;
*Save the values of gamma, lambda and the lower bound on the LR of the current best iteration
    gammaBest=gamma;
    lambdaBest=lambdaprevious;
    lbLR = bound_lr.l;
    );

* Check convergence
convergence=0;
deltalambda = abs(lambdaprevious-lambda) ;
if( deltalambda < 0.0001, contin = 0; display 'lambdas same'; convergence = 1 );

* Results output
results(iter,'deltalambda') = deltalambda;
results(iter,'noimprov') = noimprovement;
results(iter,'theta') = theta;
results(iter,'step') = stepsize;
results(iter,'gamma') = gamma ;
results(iter,'lambda') = lambda ;
results(iter,'gap') = (((-lowerbound)+upperbound)/(-lowerbound))    ; 
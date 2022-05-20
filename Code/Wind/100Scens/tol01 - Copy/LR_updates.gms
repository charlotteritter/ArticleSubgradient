*** UPDATES FOR LAMBDA, GAMMA, THETA

* Gamma update
gamma =  threshold - sum(scen, last_z(scen));

scalar it;
* Bound and theta Update -> in stepsize 6
if (bound > lowerbound,
         lowerbound = bound;
         noimprovement = 0;
         it=iter; 
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

scalar BBP;
if (it=iter,
    if (gamma ge 0 and lambda*gamma eq 0,
        BBP=1;
        );
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
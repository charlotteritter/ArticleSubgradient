*** Guarantees that never zero in the denominator
*norm = max(abs(gamma), 1e-6);
norm2 = max(gamma*gamma, 1e-6);

*** constant step-size rule
if(steprule=1,stepsize =0.015;);   

*** constant step length
if(steprule=2,stepsize = 1.5/norm2;);

*** square summable but not summable
if(steprule=3,stepsize =0.03/(3+num_iter););

*** nonsummable diminishing
if(steprule=4,stepsize= 0.005/sqrt(num_iter););

*** nonsummable diminishing step length
if(steprule=5,stepsize= (6/num_iter)/norm2;);

*** original
if(steprule=6,stepsize=theta*(upperbound - bound)/norm2;);
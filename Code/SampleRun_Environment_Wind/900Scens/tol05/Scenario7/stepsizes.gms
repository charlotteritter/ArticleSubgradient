*** Guarantees that never zero in the denominator
*norm = max(abs(gamma), 1e-6);
norm2 = max(gamma*gamma, 1e-6);

*** constant step-size rule
if(steprule=1,stepsize =0.58;);   

*** constant step length
if(steprule=2,stepsize = 200/sqrt(norm2););

*** square summable but not summable
if(steprule=3,stepsize =3/(2+num_iter););

*** nonsummable diminishing
if(steprule=4,stepsize= 1/sqrt(num_iter););

*** nonsummable diminishing step length
if(steprule=5,stepsize= (300/num_iter)/sqrt(norm2););

*** original
if(steprule=6,stepsize=theta*(upperbound - bound)/norm2;);
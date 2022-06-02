set iter /iter1*iter30/;
parameter profit(iter);
loop(iter,
    profit(iter)=1;
    );
profit('iter1')=3;
scalar a;
a=smin(iter, profit(iter));
display a;
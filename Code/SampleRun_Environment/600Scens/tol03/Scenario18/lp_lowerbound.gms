start_time = jnow ;
options optca=0 ;
options optcr=0 ;
solve schedule using RMIP minimizing OBJ ;
end_time =jnow ;
LP_time = ghour(end_time - start_time)*3600 + gminute(end_time - start_time)*60 + gsecond(end_time - start_time);
init_lambda  = Const1_2.m ;
lowerbound   = Obj.l ;
scalar LP_bound ;
lambda       = init_lambda ;
display init_lambda, LP_time, lowerbound    ;
LP_bound   = Obj.l ;

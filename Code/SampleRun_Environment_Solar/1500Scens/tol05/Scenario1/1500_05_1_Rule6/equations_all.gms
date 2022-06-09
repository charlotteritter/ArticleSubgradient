POSITIVE VARIABLES P(scen,t), Q(scen,t), Y(T), X(scen,t) ;
VARIABLES OBJ, BOUND_LR ;
BINARY VARIABLE Z(scen) ;

scalar counter ;

EQUATIONS
        Objective
        Const1(scen,t)    balance constraint
        Const2(scen,t)    max charge
        Const3(scen,t)    max discharge
        Const_chance_1(scen,t)    chance constraint big M
        Const_chance_PH(scen,t)
        Const_chance_2            chance constraint sum probabilities
        LR                lagrangian relaxation objective
        Objective_scenario(scen)
        Const1_scenario(scen,t)    balance constraint single scenario
        Const2_scenario(scen,t)    max charge    single scenario
        Const3_scenario(scen,t)    max discharge single scenario
        Const_chance_1_scenario(scen,t)    chance constraint big M  single scenario
         Const_chance_2_scenario(scen,t)
        LR_scenario(scen)
        LR_scenario_2(scen)
        LR_lb(scen)
        ;

Objective.. OBJ=E= SUM(T,Prices(T, 'REW')*Y(T) - probability*Sum(scen, Prices(T, 'CHAR')* P(scen,t) + Prices(t, 'DISCHAR') * Q(scen,t)  ) )    ;

Const1(scen,t)$(ord(t) lt card(t))..
         X(scen,t+1) =E= X(scen,t) + eta* P(scen,t) - (1/eta)* Q(scen,t) ;
Const1_scenario(scen,t)$(ord(t) lt card(t) and (ord(scen) eq counter))..
         X(scen,t+1) =E= X(scen,t) + eta* P(scen,t) - (1/eta)* Q(scen,t) ;


Const_chance_1(scen,t).. Y(T) + P(scen,t) -  Q(scen,t) -SOLAR(scen,t) =L= Z(scen)*BigM(scen,t) ;
Const_chance_1_scenario(scen,t)$(ord(scen) eq counter)..
         Y(T) + P(scen,t) -  Q(scen,t) -SOLAR(scen,t) =L=0 ;
Const_chance_2_scenario(scen,t)$(ord(scen) eq counter)..
         Y(T) + P(scen,t) -  Q(scen,t) -SOLAR(scen,t) =L= Z(scen)*BigM(scen,t) ;


Const_chance_2..      - sum(scen, z(scen)) =G= -threshold;

LR.. bound_lr =e=   SUM(T,Prices(T, 'REW')*Y(T) - probability*Sum(scen, Prices(T, 'CHAR')* P(scen,t) + Prices(t, 'DISCHAR') * Q(scen,t)  ) )
                         - lambda* (threshold - sum(scen, z(scen)))  ;

Objective_scenario(scen)$(ord(scen) eq counter)..
         OBJ =E= SUM(T,Prices(T, 'REW')*Y(T) - Prices(T, 'CHAR')* P(scen,t) + Prices(t, 'DISCHAR') * Q(scen,t)  )     ;


LR_scenario(scen)$(ord(scen) eq counter)..
         bound_lr =e=   SUM(T,Prices(T, 'REW')*Y(T) - (Prices(T, 'CHAR')* P(scen,t) + Prices(t, 'DISCHAR') * Q(scen,t))  )  + lambda* z(scen)/probability ;
LR_scenario_2(scen)$(ord(scen) eq counter)..
         bound_lr =e=   SUM(T,Prices(T, 'REW')*Y(T) - (Prices(T, 'CHAR')* P(scen,t) + Prices(t, 'DISCHAR') * Q(scen,t))  )  + lambda* z(scen)/probability +
                                  sum(t,weight_previous(scen,t)*Y(t)) + 0.5*sum(t, rho(t)*power(y_average_previous(t) - Y(t),2)) ;

LR_lb(scen)$(ord(scen) eq counter)..
         bound_lr =e=   SUM(T,Prices(T, 'REW')*Y(T) - (Prices(T, 'CHAR')* P(scen,t) + Prices(t, 'DISCHAR') * Q(scen,t))  )  + lambda* z(scen)/probability +
                                 sum(t,weight_previous(scen,t)*Y(t)) ;

*** bounds on any variables
x.up(scen,t) = BigX ;
x.lo(scen,t) = LowX ;
q.up(scen,t) = max_discharge ;
p.up(scen,t) =  max_charge ;
x.fx(scen,'t1') = X_0 ;
z.prior(scen)   = 1;


parameter last_x(scen,t), last_p(scen,t), last_q(scen,t), last_z(scen), last_ph(scen) ;
******* ALL MODELS

model schedule     / Objective,  Const1, Const_chance_1, Const_chance_2/ ;
model schedule_scenario     / Objective_scenario,  Const1_scenario,  Const_chance_1_scenario/ ; 
model INITIAL               / LR_scenario,         Const1_scenario,  Const_chance_2_scenario/ ;
model Lagrangian      / LR,    Const1, Const_chance_1/ ;

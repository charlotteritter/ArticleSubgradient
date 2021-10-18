POSITIVE VARIABLES X(W,T), Y(T);
VARIABLES OBJ, BOUND_LR;
BINARY VARIABLE Z(SCEN), R(W,T), U(W,T), V(W,T) ;

scalar counter ;

EQUATIONS
        Objective
        Const1_1(scen,t)    
        Const1_2    
        Const1_2_new
        Const_3_1(scen,t)    
        Const_3_1_scenario(scen,t)    
        Const_3_2(scen,t)
        Const_3_2_scenario(scen,t)
        Const_4_1(scen,t)    
        Const_4_1_scenario(scen,t)
        Const_4_2(scen,t)    
        Const_4_2_scenario(scen,t)
        Const_5(scen,t)
        Const_5_scenario(scen,t)
        Const_6(scen,t)
        Const_6_scenario(scen,t)
        LR          
        Objective_scenario(scen)
        Const1_1_scenario(scen,t)   
        ;

Objective.. OBJ=E= SUM(T,PRICES(T, 'REW')*Y(T) - PROBABILITY*SUM(scen, PRICES(T, 'GEN')*( X(scen,T) + GG*R(scen,t) ) )  )  ;

LR.. bound_lr =e=   SUM(T,PRICES(T, 'REW')*Y(T) - PROBABILITY*SUM(scen, PRICES(T, 'GEN')*( X(scen,T) + GG*R(scen,t) ) )  )
                         - lambda* (threshold - sum(scen, z(scen)))  ;
                         
Objective_scenario(scen)$(ord(scen) eq counter)..
         OBJ =E= SUM(TT,PRICES(tt, 'REW')*Y(TT) - PRICES(TT, 'GEN')*( X(scen,TT) + GG*R(scen,tt) ) )     ;

*Const1_1(W,T).. Y(T)-X(W,T)-R(W,T)*GG -Z(W)*WIND(W,T) =L= (1-Z(W))*BigMM(w) ;
Const1_1(scen,T).. Y(T)-X(scen,T)-R(scen,T)*GG -WIND(scen,T) =L= Z(scen)*BigM(scen,t) ;

*Const1_1_scenario(scen,t)$(ord(t) lt card(t) and (ord(scen) eq counter))..
Const1_1_scenario(scen,t)$(ord(scen) eq counter)..

         Y(T)-X(scen,T)-R(scen,T)*GG -WIND(scen,T) =L= 0 ;

Const1_2..   -SUM(SCEN, Z(SCEN)) =G= -threshold;
Const1_2_new.. SUM(SCEN, Z(SCEN)) =L= threshold;

****************************************************************************TODO Constraints (all but no. 1) in scenario form (z(w)=0)

* The generator constraints

*Ramp

Const_3_1(scen,T)$( ord(t) le (card(t)-1)).. X(scen,T+1) + GG*R(scen,t+1)- X(scen,T) - GG*R(scen,t)  =L= ramp*(U(scen,t+1) + R(scen,t)) ;
Const_3_1_scenario(scen,T)$( ord(t) le (card(t)-1) and (ord(scen) eq counter)).. X(scen,T+1) + GG*R(scen,t+1)- X(scen,T) - GG*R(scen,t)  =L= ramp*(U(scen,t+1) + R(scen,t)) ;

Const_3_2(scen,T)$( ord(t) le (card(t)-1)).. X(scen,T) + GG*R(scen,t)- X(scen,T+1) -GG*R(scen,t+1)=L= ramp*(V(scen,t+1) + R(scen,t+1)) ;
Const_3_2_scenario(scen,T)$( ord(t) le (card(t)-1) and (ord(scen) eq counter)).. X(scen,T) + GG*R(scen,t)- X(scen,T+1) -GG*R(scen,t+1)=L= ramp*(V(scen,t+1) + R(scen,t+1)) ;

*On/off

Const_4_1(scen,T)$(ord(t) ge uptime)..   Sum(tt$((ORD(TT) le ORD(T)) and (ORD(TT) ge (ORD(T) - UPTIME +1))),U(scen,TT)) =L= R(scen,T) ;
Const_4_1_scenario(scen,T)$(ord(t) ge uptime and (ord(scen) eq counter))..   Sum(tt$((ORD(TT) le ORD(T)) and (ORD(TT) ge (ORD(T) - UPTIME +1))),U(scen,TT)) =L= R(scen,T) ;

Const_4_2(scen,T)$(ord(t) ge downtime )..Sum(tt$((ORD(TT) le ORD(T)) and (ORD(TT) ge (ORD(T) - DOWNTIME +1))),V(scen,TT)) =L= 1-R(scen,T) ;
Const_4_2_scenario(scen,T)$(ord(t) ge downtime and (ord(scen) eq counter))..Sum(tt$((ORD(TT) le ORD(T)) and (ORD(TT) ge (ORD(T) - DOWNTIME +1))),V(scen,TT)) =L= 1-R(scen,T) ;

Const_5(scen,T)$(ord(t) ge 2).. U(scen,T) - V(scen,T) =E= R(scen,T) - R(scen,T-1) ;
Const_5_scenario(scen,T)$(ord(t) ge 2 and (ord(scen) eq counter)).. U(scen,T) - V(scen,T) =E= R(scen,T) - R(scen,T-1) ;

*if uptime more than 1
Const_6(scen,T)$( ord(t) le (card(t)-1))..   (G-GG)*R(scen,T)- (G-ramp)*U(scen,T) - (G-ramp)*V(scen,T+1) =G= X(scen,T);
Const_6_scenario(scen,T)$( ord(t) le (card(t)-1) and (ord(scen) eq counter))..   (G-GG)*R(scen,T)- (G-ramp)*U(scen,T) - (G-ramp)*V(scen,T+1) =G= X(scen,T);


*** bounds on any variables
X.UP(scen,T) =  G-GG;
U.UP(scen,T) = 1;
V.UP(scen,T) =1;
Z.UP(scen)=1;

* initialize the on/off variables
* assume the unit was on for last (uptime -1) periods
r.fx(scen,t)$(ord(t) eq 1) =1;
u.fx(scen,t)$(ord(t) eq 1) = 0;
v.fx(scen,t)$(ord(t) eq 1) = 0 ;
* assume the generator was producing minimim power in last time period
x.up(scen,t)$(ord(t) eq 1)= ramp - GG;
        

******* ALL MODELS

model schedule     / Objective,  Const1_1, Const1_2, Const_3_1, Const_3_2, Const_4_1, Const_4_2, Const_5, Const_6 /;
model schedule_scenario     / Objective_scenario,  Const1_1_scenario, Const_3_1_scenario, Const_3_2_scenario, Const_4_1_scenario, Const_4_2_scenario, Const_5_scenario, Const_6_scenario / ; 
*model INITIAL               / / ;
model Lagrangian      / LR,    Const1_1, Const_3_1, Const_3_2, Const_4_1, Const_4_2, Const_5, Const_6 /;

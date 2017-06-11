function f = NG(x)
%% Initial concentrations of proteins estimated from literature
AGT = 17030e-6;
ANGI = 270e-9;
ANGII = 21e-9;
ANG1_7 = 120e-9;
ANG1_9 = 60e-9;
ANGIII = 11e-9;
AT1 = 45e-9;
AT2 = 19e-9;
%% Coefficients from literature
hAGT = log(2)/(10*60);
hANGI = log(2)/(0.62/60);
hANGII = log(2)/(18/60);
hANG1_7 = log(2)/(30);
hANG1_9 = log(2)/(24);
hANGIII = log(2)/(0.5);
hAT1 = log(2)/(1.5);
hAT2 = log(2)/(1.5);
Km = 4.320e-8;
Vm = 3.023e-7;

%% Steady state mass balances for each species
f(1) = x(1)-(Vm/(Km+AGT)-hAGT)*AGT; % dAGT/dt
f(2) = Vm/(Km+AGT)*AGT-((x(2) + x(3) + x(4) + x(5) + hANGI))*ANGI; %dANGI/dt
f(3) = (x(2) + x(3))*ANGI-(x(5) + x(6) + x(7) + x(8) + hANGII)*ANGII; %dANGII/dt
f(4) = x(4)*ANGI+x(5)*ANGII-hANG1_7*ANG1_7; %dANG1_7/dt
f(5) = x(5)*ANGI-hANG1_9*ANG1_9; %dANG1_9/dt
f(6) = x(6)*ANGII-hANGIII*ANGIII; %dANGIII/dt
f(7) = x(7)*ANGII-hAT1*AT1; %dAT1/dt
f(8) = x(8)*ANGII-hAT2*AT2; %dAT2/dt 

x0 = [0,0,0,0,0,0,0,0];
x = fsolve(@NG,x0)
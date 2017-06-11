 function param_estimation
%% Take data from plots
Xdata = [1:8]; % GLU

Ydata = [1:8];% ANGII
% Ydata = [17030e-6,270e-9, 21e-9, 120e-9, 60e-9, 11e-9, 45e-9, 19e-9]
format short e
opts = optimoptions(@fsolve, 'TolFun', 1e-13,'MaxFunEvals',3000);
%'levenberg-marquardt',0.01);
% NGcoefguess = [2.11e-5,0.2655,0.07965,0.0101,0.006,0.6,1.2,0.396]; %with
% kAGT as coef(1)
%1.8138e-5
AGT = 17030e3;
AGT0 = AGT;
ANGI = 270;
ANGII = 21;
ANG1_7 = 120;
ANG1_9 = 60;
ANGIII = 11;
AT1 = 45;
AT2 = 19;
%% Coefficients for half life from literature
hAGT = log(2)/(10*3600);
hANGI = log(2)/(0.62);
hANGII = log(2)/(18);
hANG1_7 = log(2)/(30*60);
hANG1_9 = log(2)/(24*60);
hANGIII = log(2)/(0.5*60);
hAT1 = log(2)/(1.5*60);
hAT2 = log(2)/(1.5*60);
Km = 1250;
%% initial guesses
% c_ace2=hANG1_9*ANG1_9/ANGI;
% c_at1=hAT1*AT1/ANGII; %dAT1/dt
% c_at2=hAT2*AT2/ANGII;
% c_apa=hANGIII*ANGIII/ANGII;
% c_nep=(-c_ace2*ANGII+hANG1_7*ANG1_7)/ANGI;
% c_ace=(c_apa+c_ace2+c_at1+c_at2+hANGII)*ANGII/ANGI/(1+1.5/5);
% c_nonace = c_ace*1.5/5;
% Vm=((c_ace + c_nonace + c_nep + c_ace2 + hANGI))*ANGI*(AGT+Km)/AGT;
% kAGT = (Vm*AGT/(Km+AGT0)) + (hAGT*AGT);
NGcoefguess = ones(1,9);%[Vm,c_nep,c_ace2,c_apa,c_at1,c_at2,c_ace,c_nonace,kAGT] 

coef=fsolve(@NG,NGcoefguess,opts)
Vm = coef(1)
c_nep = coef(2)
c_ace2 = coef(3)
c_apa = coef(4)
c_at1 = coef(5)
c_at2 = coef(6)
c_ace = coef(7)
c_nonace = coef(8)
kAGT = coef(9)
Vm/(AGT+Km)
% Vm/(Km+AGT0)*AGT-( + c_nep + c_ace2 + hANGI)*ANGI
% kAGT-(Vm/(Km+AGT)+hAGT)*AGT

% [T,AGT] = ode45(@(t,AGT) AGTder(t,AGT,kAGT,Vm,Km,AGT0,hAGT),[0:10:10^8],17030e3);
% plot(T,AGT)
%% Set optimization options and initial guess for coefficients
% OPTIONS = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','TolX', 1e-14, 'TolFun', 1e-14); %can use optimset instead of optimoptions
% LB = -10;%1e-14; % lower bound
% UB = 600;%1e-1; % upper bound

coefficientsguess = coef; % initial guess for all coefficients
%% Parameter estimation routine lsqcurvefit
coefficients = lsqcurvefit(@glucoseRASss, coefficientsguess, Xdata, Ydata)%,LB,UB,OPTIONS); %the function returns ANGII concentration at steady state

Ycalc = glucoseRASss(coefficients,Xdata);
%% Plots
figure(2)
plot(Xdata,Ydata,'o')
hold on
plot(Xdata,Ycalc)
hold off
legend(['Data'],['Fit with coefficients'])
     function dAGTdt = AGTder(t,AGT,kAGT,Vm,Km,AGT0,hAGT)
         dAGTdt =kAGT-(Vm/(Km+AGT0)+hAGT)*AGT;
     end

end



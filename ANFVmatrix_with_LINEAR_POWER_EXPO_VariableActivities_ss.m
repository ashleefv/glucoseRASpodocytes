% Glucose(i) from 5 to 40 (mM)
% glucosevector = [5, 25,40];% 5:1:40
%  out = zeros(size(glucosevector));
% x = zeros(8,length(glucosevector));
% % PRA = Vm*[AGT]/(Km+AGT0)
AGT0 = 17030*10^-6; % mol/L
Km = 1.25E-6 ;
%  for j = 1:length(glucosevector)
%      i = glucosevector(j);
    % Units in mol/L/min
  
    kAGT = 6.3138e-7 %2.118*10^-5;
    %% variable-LINEAR
    %[fitresult,gof] = fit(x,y,model);
    Vm = 3.0351e-7 %((1.554*10^-7) * i + (1.736*10^-5));
    c_ace = 4.4387e-03 %((3.342*10^-3) * i + 0.2487); % Units in 1/min
    c_nonace = 1.3316e-03 %(0.008657*i + 0.03621);
    %% EXPONENTIAL
    % PRA = ((1.747*10^-5) * exp ((7.47*10^-3)*(i))) ;
    % c_ace = (0.252*exp(0.01043*i));
    % c_nonace = (0.06352*exp(0.04488*i));
    %% POWER LAW
    % PRA = (1.4805*10^-5)*i^(0.12617) ;
    % c_ace = (0.20014)*i^(0.17558);
    % c_nonace = (0.02326)*i^(0.7555);
    %% LOGARITHMIC
    % PRA = ((2.6168*10^-6)*log(i))+(1.3926*10^-5) ;
    % c_ace = (0.05626)*log(i)+0.175;
    % c_nonace = (0.1457*log(i))-0.15501;
    %% Other activity coefficients-constants
    c_nep = 1.6283e-04%0.0101;
    c_ace2 = 1.0697e-04 %0.006;
    c_apa = 1.2103e-02 %0.6;
    c_at1 = 1.6504e-02 %1.2;
    c_at2 =  6.9681e-03  %0.396;
    % Half life parameters=  -(ln(2)/halflife)
    % Units in 1/sec
    hAGT = log(2)/(10*3600);
    hANGI = log(2)/(0.62);
    hANGII = log(2)/(18);
    hANG1_7 = log(2)/(30*60);
    hANG1_9 = log(2)/(24*60);
    hANGIII = log(2)/(0.5*60);
    hAT1 = log(2)/(1.5*60);
    hAT2 = log(2)/(1.5*60);
    %A*x=b
    % A = [coefficients]
    % x = [concentrations]
    % b = [rate]
    A = [-(Vm/(Km+AGT0))-hAGT,0,0,0,0,0,0,0; ...
        Vm/(Km+AGT0),-(c_ace + c_nonace + c_nep + c_ace2 + hANGI),0,0,0,0,0,0;...
        0,(c_ace + c_nonace),-(c_ace2 + c_apa + c_at1 + c_at2 + hANGII),0,0,0,0,0;...
        0,c_nep,c_ace2,(-hANG1_7),0,0,0,0; ...
        0,c_ace2,0,0,(-hANG1_9),0,0,0;...
        0,0,c_apa,0,0,-hANGIII,0,0;...
        0,0,c_at1,0,0,0,-hAT1,0;...
        0,0,c_at2,0,0,0,0,-hAT2];
    b = [(-kAGT); 0; 0; 0; 0; 0; 0; 0];
    format shortEng
    format compact 
%      x(:,j) = A\b
 x = A\b

% glucosevector
% 5out = x(3,:)./x(3,1)
% 
 % plot(glucosevector,x(6,:),'o')
 % hold on
 % plot(glucosevector,x(3,:),'x')
 % hold off


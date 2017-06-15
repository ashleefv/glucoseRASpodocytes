function out = glucoseRASssA3(coefficients,GLU,baseline,scenario,c_nep,c_ace2,c_apa,c_at2,printoutput)

for i = 1:length(GLU)
    % The scenarios intpret glucose dependent functional forms with
    % variable numbers of input coefficients to transform to the 9 needed
    % coefficients to solve the glucoseRASss model

 coef_output = Approach3(coefficients, GLU(i));

kAGT = coef_output(1);
c_ACE = coef_output(2);
c_at1 = coef_output(3);
Vm = coef_output(4);
%% Coefficients from literature
S = 17030e3;
Km = 1250;
hAGT = log(2)/(10*3600);
hANGI = log(2)/(0.62);
hANGII = log(2)/(18);
hANG1_7 = log(2)/(30*60);
hANG1_9 = log(2)/(24*60);
hANGIII = log(2)/(0.5*60);
hAT1 = log(2)/(1.5*60);
hAT2 = log(2)/(1.5*60);
% glucoseRASss is a function describing glucose-sensitive RAS in podocytes
% at steady state
% the function returns ANGII concentration at steady state
% the input to the function is the Xdata value for glucose concentration
% and the current values for the changing coefficients
A = [-hAGT-Vm,0,0,0,0,0,0,0; ...
        Vm,-(c_ACE + c_nep + c_ace2 + hANGI),0,0,0,0,0,0;...
        0,(c_ACE),-(c_ace2 + c_apa + c_at1 + c_at2 + hANGII),0,0,0,0,0;...
        0,c_nep,c_ace2,(-hANG1_7),0,0,0,0; ...
        0,c_ace2,0,0,(-hANG1_9),0,0,0;...
        0,0,c_apa,0,0,-hANGIII,0,0;...
        0,0,c_at1,0,0,0,-hAT1,0;...
        0,0,c_at2,0,0,0,0,-hAT2];

b = [(-kAGT); 0; 0; 0; 0; 0; 0; 0];
%options for coefficents
x = A\b;
%linear algebra matrices
%define A matrix and B vector here
% coefficientsguess = [
%X=A\B;

% <<<<<<< HEAD

 x(3) ;
out(i) = x(3)/baseline; % ANGII concentration at steady state
% =======
% >>>>>>> f80245419646f6d78dab3c239234f3a9555f6ac2
if printoutput == 1
    scenario;
    coef_output;
end
if printoutput == 2
    for j = 1:8
        % all output concentration at steady state for plotting these
        % concentrations
        out(i,j) = x(j); 
    end
else
    % ANGII concentration at steady state for parameter estimation
    out(i) = x(3)/baseline; 
end
end
end
function param_estimation_Approach3
close all
% Comment out close all if you want 
%to overlapp the plots for peptide concentrations from
% fmincon code with the output from this code
format long e
%% Take data from plots
Xdata = [5, 25, 40]; % GLU %mmol/L
Ydata = [1, 1.4, 2.1];% ANGII
e = [0.09, 0.08, 0.58];
%figure(1)
l = 1 %to generate ANG II concentration plot in figure 11
h = errorbar(Xdata,Ydata,e,'kx','Markersize',12,'Linewidth',2);
s = h.LineWidth;
h.LineWidth = 1.5;
%% initial guesses
NG;
load('NGvalues','NGvalues');

c_nep = NGvalues(4);
c_ace2 = NGvalues(5);
c_apa = NGvalues(6);
c_at2 = NGvalues(8);
% if you want to print the 9 coefficents for glucoseRASssA3, printoutput = 1
% else printoutput = 0
% if you want the total output from glucoseRASssA3 output as a vector,
% printoutput = 2
% leave as default 0. Specified otherwise as needed for plotting or
% printing to command window in the code.
printoutput = 0;
baseline = glucoseRASssA3(NGvalues,5,1,0,c_nep,c_ace2,c_apa,c_at2,printoutput);

options = optimoptions('fmincon','Algorithm','sqp','TolFun',1e-16,'TolX',1e-16,'MaxFunEvals',20000,'StepTolerance',1e-16);
%% yy = 5 is Approach 3
yy = 5
scenario = yy;
%% Parameter estimation routine lsqcurvefit

% %Vm= Vma*GLU + Vmb; c_ace=c_aceA*GLU + c_aceB;c_at1=c_at1A*GLU + c_at1B
      coefficientsguess= [ 1.177568154554772e-01     1.146519658757140e-01,...
          2e-12        3.033621116991067e+02,...
        3.432927517238616e-03     3.345022358146005e-03                         2e-12,...
          3.784671163979920e+00];
     LB = zeros(size(coefficientsguess));


%% Parameter estimation routine lsqcurvefit

UB = [];

    fun = @(coefficients)sum((Ydata - glucoseRASssA3(coefficients,Xdata,baseline,scenario,c_nep,c_ace2,c_apa,c_at2,printoutput)).^2);
    A = [];
    b = [];
    Aeq = [];
    beq = [];
 %lb = [];
% ub = [];
    lb = zeros(size(coefficientsguess));
    ub = Inf(size(coefficientsguess));
    coefficients = fmincon(fun,coefficientsguess,A,b,Aeq,beq,lb,ub,[],options)



X = Xdata(1):0.5:Xdata(end);
% Y = zeros(size(X));
% Ycalc = zeros(size(Xdata));
for i = 1:length(X)
    Y(i)= glucoseRASssA3(coefficients,X(i),baseline,scenario,c_nep,c_ace2,c_apa,c_at2,printoutput);
end
for i = 1:length(Xdata)-1
    Ycalc(i)= glucoseRASssA3(coefficients,Xdata(i),baseline,scenario,c_nep,c_ace2,c_apa,c_at2,printoutput);
end
Ycalc(i+1)= glucoseRASssA3(coefficients,Xdata(end),baseline,scenario,c_nep,c_ace2,c_apa,c_at2,1);

% All output from glucoseRASssA3
for i = 1:length(X)
    Y_all(i,:)= glucoseRASssA3(coefficients,X(i),1,scenario,c_nep,c_ace2,c_apa,c_at2,2);
end
%% Plots
figure(11)
box on
get(gca);set(gca,'FontSize',10)%,'FontWeight','Bold');
ax = gca;
ax.YColor = 'black';
 plot(X(1:3:end),Y_all((1:3:end),3),'o-r','LineWidth',2,'MarkerSize',5)
xlabel('Glucose (mM)','FontSize',10)%,'FontWeight','Bold')
ylabel('ANGII (nM)','FontSize',10)

legend('Approach 3', 'Location','Southeast')



 

 %legend('Data','Approach 3','Location','Northwest')%,'Color',[0 0.5 0],'LineWidth',2,'MarkerSize',12)



%% RMSE and Norm Calculation 
RMSE(yy+1) = sqrt(mean((Ydata - Ycalc).^2))

Norm_calc(yy+1) = sum((Ydata - Ycalc).^2);

figure(1)
box on
get(gca);set(gca,'Fontsize',10);%,'FontWeight','Bold');
hold on
axis([5 40 0 4])
plot(X(1:3:end),Y(1:3:end),'o-r','LineWidth',2,'MarkerSize',5)
 legend('Data from Durvasula and Shankland, 2008',['Approach 3: RMSE = ' num2str(round(RMSE(5+1),8))], 'Location','NorthWest')
xlabel('Glucose (mM)','Fontsize',10)%,'FontWeight','Bold')
ylabel('ANGII Change from Baseline','Fontsize',10)%,'FontWeight','Bold')
hold off
axis([5 40 0.8 3.6])
 

end






   

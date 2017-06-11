function param_estimation_scenario_0_5
%close all
% Comment out close all if you want 
%to overlapp the plots for peptide concentrations from
% fmincon code with the output from this code
format long e
%% Take data from plots
Xdata = [5, 25, 40]; % GLU %mmol/L
Ydata = [1, 1.4, 2.1];% ANGII
e = [0.09, 0.08, 0.58];
% figure(1)

l = 1
%plot(Xdata,Ydata,'kx','Markersize',12)

h = errorbar(Xdata,Ydata,e,'kx','Markersize',12)
s = h.LineWidth;
h.LineWidth = 1.5
%% initial guesses
NG;
load('NGvalues','NGvalues');

c_nep = NGvalues(4);
c_ace2 = NGvalues(5);
c_apa = NGvalues(6);
c_at2 = NGvalues(8);
% if you want to print the 9 coefficents for glucoseRASss05, printoutput = 1
% else printoutput = 0
% if you want the total output from glucoseRASss05 output as a vector,
% printoutput = 2
% leave as default 0. Specified otherwise as needed for plotting or
% printing to command window in the code.
printoutput = 0;
baseline = glucoseRASss05(NGvalues,5,1,0,c_nep,c_ace2,c_apa,c_at2,printoutput);

options = optimoptions('fmincon','Algorithm','sqp','TolFun',1e-16,'TolX',1e-16,'MaxFunEvals',20000,'StepTolerance',1e-16);
for yy = 4
    
scenario = yy; %edit for desired scenario

%% Parameter estimation routine lsqcurvefit
if scenario == 1 %Vm= Vma*GLU + Vmb
    coefficientsguess = [2.805953543653534e-03    -2.757657249946220e-04 ,...
     -2.348755561601649e-02  ,...
    -1.668986022883667e+02    -1.685774047241159e+00];
%     coefficientsguess = [1.103822324074590e+00     1.101016998005601e+00 ,...
%      2e-12  ,...
%     4.950685540241108e-01     1.961355230418327e-02]; % initial guess for all coefficients
    LB = zeros(size(coefficientsguess));
elseif scenario == 2 %c_ace=c_aceA*GLU + c_aceB
    %coefficientsguess = ones(1:5);
    coefficientsguess = [   -3.609928272726562e-03  -1.240811257988384e-04  -6.715214044691388e-03     ,...
       -3.767179069478035e+00    ,...
  3.034683893640923e+02    ];
    LB = zeros(size(coefficientsguess));
elseif scenario == 3 %c_at1=c_at1A*GLU + c_at1B
    coefficientsguess =  [       3.066589908448575e-01     3.035519149069141e-01    ,...
          8.931547418727733e-01     ,...
        -5.073423330714726e-02  3.035071407475688e+02 ];
    LB = zeros(size(coefficientsguess));
elseif scenario == 4 %c_nonace=c_nonaceA*GLU + c_nonaceB
    coefficientsguess= [ 6.313823632053241e+02     4.533794480918563e-03    ,...
       1.296703909210782e-02,...
       1.705688364046031e-05   2.472978807773762e-04,...
    7.072930413876994e-04     1.527482117056147e-07]
%      coefficientsguess= [ 6.313823632053239e+02     2.057956338843580e+00     ,...
%     4.504747645041913e+00        1.705688364046031e-05,...
%      3.615305143993680e+00         -2.360527261714031e-01,...
%     -1.860363491152925e+00];
     LB = zeros(size(coefficientsguess));
elseif scenario == 5
 %coefficientsguess = ones(1:8);
 %[NGvalues(1:9),1,1,1,1];
      coefficientsguess= [ 1.177568154554772e-01     1.146519658757140e-01,...
          2e-12        3.033621116991067e+02,...
        3.432927517238616e-03     3.345022358146005e-03                         2e-12,...
          3.784671163979920e+00];
     LB = zeros(size(coefficientsguess));
else 
    coefficientsguess = NGvalues(1:9);
    LB = zeros(size(coefficientsguess));
end
%% Parameter estimation routine lsqcurvefit
% glucoseRASss0505(coefficientsguess,Xdata(1))
% options = optimoptions('fmincon','Algorithm','trust-region-reflective')
UB = [];
if scenario == 6
    coefficients = NGHG(NGvalues);
else
    fun = @(coefficients)sum((Ydata - glucoseRASss05(coefficients,Xdata,baseline,scenario,c_nep,c_ace2,c_apa,c_at2,printoutput)).^2);
    A = [];
    b = [];
    Aeq = [];
    beq = [];
 %lb = [];
% ub = [];
    lb = zeros(size(coefficientsguess));
    ub = Inf(size(coefficientsguess));
    coefficients = fmincon(fun,coefficientsguess,A,b,Aeq,beq,lb,ub,[],options)

end

X = Xdata(1):0.5:Xdata(end);
% Y = zeros(size(X));
% Ycalc = zeros(size(Xdata));
for i = 1:length(X)
    Y(i)= glucoseRASss05(coefficients,X(i),baseline,scenario,c_nep,c_ace2,c_apa,c_at2,printoutput);
end
for i = 1:length(Xdata)-1
    Ycalc(i)= glucoseRASss05(coefficients,Xdata(i),baseline,scenario,c_nep,c_ace2,c_apa,c_at2,printoutput);
end
Ycalc(i+1)= glucoseRASss05(coefficients,Xdata(end),baseline,scenario,c_nep,c_ace2,c_apa,c_at2,1);

% All output from glucoseRASss05
for i = 1:length(X)
    Y_all(i,:)= glucoseRASss05(coefficients,X(i),1,scenario,c_nep,c_ace2,c_apa,c_at2,2);
end
%% Plots
figure(11)
box on
get(gca);set(gca,'FontSize',10)%,'FontWeight','Bold');
ax = gca;
ax.YColor = 'black';
 hold on
% axis([5 40 15 45])
if scenario == 0
 plot(X(1:5:end),Y_all((1:5:end),3),'square')
elseif scenario == 1
 plot(X(1:4:end),Y_all((1:4:end),3),'d')
 elseif scenario == 2
 plot(X(1:5:end),Y_all((1:5:end),3),'*')
 elseif scenario == 3
 plot(X(1:3:end),Y_all((1:3:end),3),'X')
 elseif scenario == 4
  plot(X(1:5:end),Y_all((1:5:end),7),'o-r','LineWidth',2,'MarkerSize',5)
  axis([5 40 0 300]);
  legend('A1','A2:S4','A3','Fontsize',10, 'Location','Northwest')
xlabel('Glucose (mM)','FontSize',10)%,'FontWeight','Bold')
ylabel('AT1R-ANGII (nM)','FontSize',10)
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 2.5 2.5]);
export_fig('C:/Users/Minu/Desktop/minu_pilvankar_CHE5110/matlab/Figures/PaperVersion/fig7AT12525', '-pdf', '-png', '-eps', '-tiff'); 
 elseif scenario == 5
%yyaxis right
 plot(X(1:3:end),Y_all((1:3:end),3),'o-r','LineWidth',2,'MarkerSize',5)
%axis([5 40 0 inf])

% axis([5 40 0 inf])
legend('Approach 3','Fontsize',10, 'Location','Southeast')
xlabel('Glucose (mM)','FontSize',10)%,'FontWeight','Bold')
ylabel('ANGII (nM)','FontSize',10)%,'FontWeight','Bold')
%legend('Approach 3','Approach 2')
hold off
 end



% 
 %legend('Scenario 0', 'Scenario 1','Scenario 2','Scenario 3', 'Scenario 4','Scenario 5','Location','Southeast')
% 
% legend('Approach 3', 'Location','Southeast')
%  set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 4 4]);
% box on
% ax = gca;
% ax.YColor = 'black';
% %get(gca);
% set(gca,'FontSize',18,'FontName','Arial');
% hold on
% % if scenario == 2
% %yyaxis left
% %plot(X(1:2:end),Y_all((1:2:end),4),'square-','LineWidth',2,'MarkerSize',14,'Color',[ 0.4660  ,  0.6740 ,   0.1880])
% plot(X,Y_all(:,3),'squareb')%axis([5 40 100 200])
% axis([5 40 15 45])
%  xlabel('Glucose (mM)','FontSize',18)
%  ylabel('ANGII (nM)','FontSize',18)
% %legend('Approach 3')
% hold off
% 
%  
% % 
%  %legend('Scenario 0', 'Scenario 1','Scenario 2','Scenario 3', 'Scenario 4','Scenario 5','Location','Southeast')
% % 
% %  legend('Approach 3', 'Location','Southeast')
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 2.5 2.5]);
%%export_fig('C:/Users/Minu/Desktop/minu_pilvankar_CHE5110/matlab/Figures/PaperVersion/fig5cA3ANGIIconc2525', '-pdf', '-png', '-eps', '-tiff'); 
% % 
% 
figure(1)
box on
get(gca);set(gca,'Fontsize',10)%,'FontWeight','Bold');
hold on
if scenario == 0
 plot(X(1:5:end),Y(1:5:end),'square','Linewidth',2)
elseif scenario == 1
 plot(X(1:4:end),Y(1:4:end),'d','Linewidth',2)
 elseif scenario == 2
 plot(X(1:5:end),Y(1:5:end),'*','Linewidth',2)
 elseif scenario == 3
 plot(X(1:3:end),Y(1:3:end),'X','Linewidth',2)
%  elseif scenario == 4
%  plot(X(1:4:end),Y(1:4:end),'+','Linewidth',2)
%  elseif scenario == 6
% xlabel('Glucose (mM)','Fontsize',18)%,'FontWeight','Bold')
% ylabel('ANGII Change from Baseline','Fontsize',18)%,'FontWeight','Bold')
% 
%  plot(X(1:3:end),Y(1:3:end),'o','LineWidth',2,'MarkerSize',10)%,'Color',[ 0.9290 ,   0.6940  ,  0.1250])
% axis([5 40 0.8 2.8])
  elseif scenario == 5
      xlabel('Glucose (mM)','Fontsize',10)%,'FontWeight','Bold')
ylabel('ANGII Change from Baseline','Fontsize',10)%,'FontWeight','Bold')
      axis([5 40 0.8 3.6])
 plot(X(1:3:end),Y(1:3:end),'o-r','LineWidth',2,'MarkerSize',5)
 legend('Data','Approach 3','Location','Northwest')%,'Color',[0 0.5 0],'LineWidth',2,'MarkerSize',12)
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 2.5 2.5]);
 %%%export_fig('C:/Users/Minu/Desktop/minu_pilvankar_CHE5110/matlab/Figures/PaperVersion/fig6cA3ANGII2525', '-pdf', '-png', '-eps', '-tiff'); 
% axis([5 40 0.8 2.8])
% xlabel('Glucose (mM)','Fontsize',18)%,'FontWeight','Bold')
% ylabel('ANGII (nM)','Fontsize',18)%,'FontWeight','Bold')

hold off
 end
hold off
RMSE(yy+1) = sqrt(mean((Ydata - Ycalc).^2))
%legend(['Approach 3: RMSE = ' num2str(round(RMSE(5+1),8))], 'Location','EastOutside')
Norm_calc(yy+1) = sum((Ydata - Ycalc).^2);
if scenario == 1 
figure(2)
subplot(221)
plot(X, coefficients(5)*X+coefficients(4))
xlabel('Glucose (mM)')
ylabel('V_m (mM/s)')
titletext = 'Scenario %i: Vm = %1.2e [Glucose] + %1.2e\n';
g= sprintf(titletext,scenario, coefficients(5),coefficients(4));
title(g)
elseif scenario == 2
figure(2)
subplot(222)
plot(X, coefficients(2)*X+coefficients(1))
xlabel('Glucose (mM)')
ylabel('c_{ace} (mM/s)')
titletext = 'Scenario %i: c_{ace} = %1.2e [Glucose] + %1.2e\n';
g= sprintf(titletext,scenario, coefficients(2),coefficients(1));
title(g)
elseif scenario == 3
figure(2)
subplot(223)
plot(X, coefficients(4)*X+coefficients(3))
xlabel('Glucose (mM)')
ylabel('c_{at1} (mM/s)')
titletext = 'Scenario %i: c_{at1} = %1.2e [Glucose] + %1.2e\n';
g= sprintf(titletext,scenario, coefficients(4),coefficients(3));
title(g)
% title(['Scenario 3: c_{at1} = ' num2str(coefficients(7)) 'Glucose + ' num2str(coefficients(9))])
elseif scenario == 4
figure(2)
subplot(224)
plot(X, coefficients(3)*X+coefficients(2))
xlabel('Glucose (mM)')
ylabel('c_{nonace} (nM/s)')
titletext = 'Scenario %i: c_{nonace} = %1.2e [Glucose] + %1.2e\n';
g= sprintf(titletext,scenario, coefficients(3),coefficients(2));
title(g)
elseif scenario == 5
figure(3)
subplot(221)
plot(X, coefficients(5)*ones(size(X)))
xlabel('Glucose (mM)')
% ylabel('Vm (mM/s)')
elseif scenario == 6
figure(3)
subplot(221)
plot(X, coefficients(10)*ones(size(X)))
xlabel('Glucose (mM)')
% ylabel('Vm (mM/s)')
end
end
% figure(1)
% axis([5 40 0.8 2.8])


% legend('Data from [1]',['Approach 3: RMSE = ' num2str(round(RMSE(5+1),8))], 'Location','NorthWest')
% legend('Data from [1]',['Case 0: RMSE = ' num2str(round(RMSE(0+1),2))],...
%     ['Case 1: RMSE = ' num2str(round(RMSE(1+1),2))],...
%     ['Case 2: RMSE = ' num2str(round(RMSE(2+1),2))],...
%     ['Case 3: RMSE = ' num2str(round(RMSE(3+1),3))],...
%     ['Case 4: RMSE = ' num2str(round(RMSE(4+1),2))],...
% figure(1)

   % ['Case 6: RMSE = ' num2str(round(RMSE(6+1),8))],


% set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 2.5 2.5]);
% %%%export_fig('C:/KidneyDiabetes/matlab/Figures/PaperVersion/A3Modelvalid2525', '-pdf', '-png', '-eps', '-tiff'); 






   

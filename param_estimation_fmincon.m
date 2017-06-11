function param_estimation_fmincon
close all
format long e
%% Take data from plots
Xdata = [5, 25, 40]; % GLU %mmol/L
Ydata = [1, 1.4, 2.1];% ANGII
e = [0.09, 0.08, 0.58];
%figure(1)
l = 7
h = errorbar(Xdata,Ydata,e,'kx','Markersize',12,'Linewidth',2);
s = h.LineWidth;
h.LineWidth = 1.5;
%% initial guesses
% opts = optimoptions(@fsolve, 'TolFun', 1e-16,'MaxFunEvals',30000);
% NGcoefguess = ones(1,8); %the matrix may change in size need to find a way to standardize
NG;
load('NGvalues','NGvalues')
%kAGT = NGvalues(2);
% if you want to print the 9 coefficents for glucoseRASss, printoutput = 1
% else printoutput = 0
% if you want the total output from glucoseRASss output as a vector,
% printoutput = 2
% leave as default 0. Specified otherwise as neede;d for plotting or
% printing to command window in the code.
printoutput = 0;
baseline = glucoseRASss(NGvalues,5,1,0,printoutput);

options = optimoptions('fmincon','Algorithm','sqp','TolFun',1e-16,'TolX',1e-16,'MaxFunEvals',95000);
for yy =  [5,4]
    
scenario = yy; %edit for desired scenario6

%% Parameter estimation routine lsqcurvefit
if scenario == 1 %Vm= Vma*GLU + Vmb
    coefficientsguess = [6.313823632053241e+02     2.805953543653534e-03     1.454338200073524e-03,...
    -3.871170456138443e-02    -2.788701138364793e-02    -2.348755561601649e-02    -3.305016330914784e-02,...
    -1.668986022883667e+02    -1.685774047241159e+00]; % initial guess for all coefficients
    LB = zeros(1,9);
elseif scenario == 2 %c_ace=c_aceA*GLU + c_aceB
    coefficientsguess = [ 6.313823632053241e+02    -3.609928272726562e-03         1.346875151901978e+00,...
    -2.441788832026719e+00    -3.771580307207119e+00    -3.767179069478035e+00    -3.776716636512934e+00,...
  3.034683893640923e+02    -1.240811257988384e-04];
    LB = zeros(1,9);
elseif scenario == 3 %c_at1=c_at1A*GLU + c_at1B
    coefficientsguess = [  6.313823632053241e+02     3.066589908448575e-01         2.884945324408734e-02,...
       9.054443071354151e-01     8.887538073930262e-01     8.931547418727733e-01     8.836193838333236e-01,...
     3.035071407475688e+02    -5.073423330714726e-02];
    LB = zeros(1,9);

elseif scenario == 4
%      coefficientsguess = [NGvalues(1:9),1,1,1,1];
coefficientsguess= [ 6.313823632053241e+02     4.533794480918563e-03     1.628277841850352e-04 ,...
    1.069671574938187e-04    1.210256981930063e-02     1.296703909210782e-02,...
    6.968146259597334e-03     1.705688364046031e-05   2.472978807773762e-04,...
    7.072930413876994e-04     1.527482117056147e-07]
% coefficientsguess= [6.313823632053239e+02     2.057956338843580e+00     ,...
%     -2.701399327671916e+00     1.786788271877614e+00     4.500347617695135e+00,...
%   4.504747645041913e+00     4.495213212605229e+00     1.705688364046031e-05,...
%      3.615305143993680e+00         -2.360527261714031e-01,...
%     -1.860363491152925e+00];
     LB = zeros(1,11);
else 
    coefficientsguess = NGvalues(1:8);
    LB = zeros(1,8);
end
%% Parameter estimation routine lsqcurvefit
% glucoseRASss(coefficientsguess,Xdata(1))
% options = optimoptions('fmincon','Algorithm','trust-region-reflective')
UB = [];
if scenario == 5
    coefficients = NGHG(NGvalues)
else
    fun = @(coefficients)sum((Ydata - glucoseRASss(coefficients,Xdata,baseline,scenario,printoutput)).^2);
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

for i = 1:length(X)
    Y(i)= glucoseRASss(coefficients,X(i),baseline,scenario,printoutput);
end
for i = 1:length(Xdata)-1
    Ycalc(i)= glucoseRASss(coefficients,Xdata(i),baseline,scenario,printoutput);
end
Ycalc(i+1)= glucoseRASss(coefficients,Xdata(end),baseline,scenario,1);


% All output from glucoseRASss
for i = 1:length(X)
    Y_all(i,:)= glucoseRASss(coefficients,X(i),1,scenario,2);
end
%% Plots
figure(11)
box on
get(gca);set(gca,'FontSize',10);
%,'FontWeight','Bold');
ax = gca;
ax.YColor = 'black';
hold on
if scenario == 0
 plot(X(1:5:end),Y_all((1:5:end),3),'square','Linewidth',2)
elseif scenario == 1
 plot(X(1:4:end),Y_all((1:4:end),3),'d','Linewidth',2)
 elseif scenario == 2
 plot(X(1:5:end),Y_all((1:5:end),3),'*','Linewidth',2)
 elseif scenario == 3
 plot(X(1:3:end),Y_all((1:3:end),3),'X','Linewidth',2)
  elseif scenario == 4
%yyaxis right
 plot(X(1:3:end),Y_all((1:3:end),l),'d-','Color',[ 0.4660    0.6740    0.1880],'LineWidth',2)
 xlabel('Glucose (mM)','Fontsize',10)%,'FontWeight','Bold')
ylabel('ANGII (nM)','Fontsize',10)
axis([5 40 0 Inf])%,'FontWe
 legend('Scenario 0',...
    'Scenario 1',...
    'Scenario 2',...
    'Scenario 3',...
    'Scenario 4','Location','Southeast')
    %set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 2.5 2.5]);
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 2.5 2.5]);
%%export_fig('C:/Users/Minu/Desktop/minu_pilvankar_CHE5110/matlab/Figures/PaperVersion/fig5bA2ANGIIconc2525', '-pdf', '-png', '-eps', '-tiff');%,'LineWidth',2,'MarkerSize',8)%,'Color',[ 0 ,   0.5  ,  0])%,'MarkerSize',8)
%axis([5 40 0 inf])

elseif scenario == 5
%yyaxis left
plot(X(1:3:end),Y_all((1:3:end),l),'x-b','Linewidth',2,'MarkerSize',5)
xlabel('Glucose (mM)','Fontsize',10)%,'FontWeight','Bold')
ylabel('ANGII (nM)','Fontsize',10)%,'FontWeight','Bold')
legend('Approach 1','Fontsize',10, 'Location','Southeast')
axis([5 40 0 Inf])

% axis([5 40 0 inf])
% xlabel('Glucose (mM)','FontSize',18)%,'FontWeight','Bold')
% ylabel('ANGII(nM)','FontSize',18)%,'FontWeight','Bold')
%legend('Approach 1','Approach 2')
hold off
 end


% 
 %legend('Scenario 0', 'Scenario 1','Scenario 2','Scenario 3', 'Scenario 4','Scenario 5','Location','Southeast')
% 
% legend('Approach 1', 'Location','Southeast')
%  set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 3 3]);
% %export_fig('C:/KidneyDiabetes/matlab/Figures/PaperVersion/A2ANGII2525','-pdf', '-png', '-eps', '-tiff'); 

figure(1)
box on
get(gca);set(gca,'Fontsize',10);%,'FontWeight','Bold');
hold on
axis([5 40 0 4])
if scenario == 0
 plot(X(1:3:end),Y(1:3:end),'square','Linewidth',2)
elseif scenario == 1
 plot(X(1:4:end),Y(1:4:end),'d','Linewidth',2)
 elseif scenario == 2
 plot(X(1:4:end),Y(1:4:end),'*','Linewidth',2)
 elseif scenario == 3
 plot(X(1:5:end),Y(1:5:end),'X','Linewidth',2)
  elseif scenario == 4
xlabel('Glucose (mM)','Fontsize',18)%,'FontWeight','Bold')
ylabel('ANGII Change from Baseline','Fontsize',18)%,'FontWeight','Bold')
axis([5 40 0.8 3.6])
 plot(X(1:5:end),Y(1:5:end),'o','Linewidth',2)%,'LineWidth',2,'MarkerSize',12)%,'Color',[ 0.9290 ,   0.6940  ,  0.1250])

  elseif scenario == 5
 plot(X(1:3:end),Y(1:3:end),'o-b','LineWidth',2,'MarkerSize',5)%,'Color',[0 0.5 0],'LineWidth',2,'MarkerSize',12)
legend('Data','Approach 1', 'Location','Northwest');
axis([5 40 0.8 3.6])
xlabel('Glucose (mM)','Fontsize',10)%,'FontWeight','Bold')
ylabel('ANGII Change from Baseline','Fontsize',10)%,'FontWeight','Bold')
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 2.5 2.5]);
%%export_fig('C:/Users/Minu/Desktop/minu_pilvankar_CHE5110/matlab/Figures/PaperVersion/fig6aA1ANGII2525', '-pdf', '-png', '-eps', '-tiff');
hold off
 end


% get(gca);set(gca,'FontSize',10,'FontName','Arial');


% set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 4 4]);
% 
%%export_fig figure(1) -r1000 -a4 -q101 -painters -eps -png -tiff
% ylabel({'ANGII Change from Baseline'},'FontSize',10)
% hold off
RMSE(yy+1) = sqrt(mean((Ydata - Ycalc).^2))

Norm_calc(yy+1) = sum((Ydata - Ycalc).^2);
if scenario == 1
figure(2)
subplot(221)
plot(X, coefficients(9)*X+coefficients(8))
xlabel('Glucose (mM)')
ylabel('V_m (mM/s)')
titletext = 'Scenario %i: Vm = %1.2e [Glucose] + %1.2e\n';
g= sprintf(titletext,scenario, coefficients(9),coefficients(8));
title(g)
elseif scenario == 2
figure(2)
subplot(222)
plot(X, coefficients(9)*X+coefficients(2))
xlabel('Glucose (mM)')
ylabel('c_{ace} (mM/s)')
titletext = 'Scenario %i: c_{ace} = %1.2e [Glucose] + %1.2e\n';
g= sprintf(titletext,scenario, coefficients(9),coefficients(2));
title(g)
elseif scenario == 3
figure(2)
subplot(223)
plot(X, coefficients(9)*X+coefficients(6))
xlabel('Glucose (mM)')
ylabel('c_{at1} (mM/s)')
titletext = 'Scenario %i: c_{at1} = %1.2e [Glucose] + %1.2e\n';
g= sprintf(titletext,scenario, coefficients(9),coefficients(6));
title(g)
% title(['Scenario 3: c_{at1} = ' num2str(coefficients(7)) 'Glucose + ' num2str(coefficients(9))])
elseif scenario == 4
figure(3)
subplot(221)
plot(X, coefficients(9)*ones(size(X)))
xlabel('Glucose (mM)')
% ylabel('Vm (mM/s)')
elseif scenario == 5
figure(3)
subplot(221)
plot(X, coefficients(10)*ones(size(X)))
xlabel('Glucose (mM)')
% ylabel('Vm (mM/s)')
end
end

% figure(1)
% axis([5 40 0.8 2.8])
% xlabel('Glucose (mM)','Fontsize',10)
 figure(1)
% %uncomment next 6 lines to give legends for scenario 0 to 5
axis([5 40 0.8 3.6])
xlabel('Glucose (mM)','Fontsize',10)%,'FontWeight','Bold')
ylabel('ANGII Change from Baseline','Fontsize',10)%,'FontWe
 legend('Data','Scenario 0',...
    'Scenario 1',...
    'Scenario 2',...
    'Scenario 3',...
    'Scenario 4','Location','Northwest')
    set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 2.5 2.5]);
 %%export_fig('C:/Users/Minu/Desktop/minu_pilvankar_CHE5110/matlab/Figures/PaperVersion/fig6bA2ANGII2525', '-pdf', '-png', '-eps', '-tiff');
   % ['Case 5: RMSE = ' num2str(round(RMSE(5+1),8))],'Location','Northwest')
%figure(1)
% legend('Data from [1]',['Scenario 0: RMSE = ' num2str(round(RMSE(0+1),2))],...
%     ['Scenario 1: RMSE = ' num2str(round(RMSE(1+1),2))],...
%     ['Scenario 2: RMSE = ' num2str(round(RMSE(2+1),2))],...
%     ['Scenario 3: RMSE = ' num2str(round(RMSE(3+1),3))],...
%     ['Scenario 4: RMSE = ' num2str(round(RMSE(4+1),2))],...
%     ['Scenario 5: RMSE = ' num2str(round(RMSE(5+1),8))],'Location','Northwest')
%legend('Data from [1]',['Approach 2: RMSE = ' num2str(round(RMSE(6+1),8))],'Location','Northwest')
      % ['Approach 2: RMSE = ' num2str(round(RMSE(5+1),8))],'Location','Northwest')
%% Uncomment the set and %export_fig command to %export the plot.
% Make sure the file name and folder is correct
% set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 2.5 2.5]);
% %export_fig('C:/KidneyDiabetes/matlab/Figures/PaperVersion/A1Modelvalid2525', '-pdf', '-png', '-eps', '-tiff'); 
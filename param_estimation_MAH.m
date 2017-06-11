function param_estimation_MAH
format long e
%% Take data from plots
Xdata = [5, 25, 40]; % GLU %mmol/L
Ydata = [1, 1.4, 2.1];% ANGII
figure(1)
% plot(Xdata,Ydata,'x','Markersize',12)
%% initial guesses
opts = optimoptions(@fsolve, 'TolFun', 1e-16,'MaxFunEvals',3000);
NGcoefguess = ones(1,9); %the matrix may change in size need to find a way to standardize
[NGvalues,FVAL,EXITFLAG,OUTPUT]=fsolve(@NG,NGcoefguess,opts);
kAGT = NGvalues(1);
printoutput = 0;
baseline = glucoseRASss(NGvalues,5,1,0,kAGT,printoutput);

for yy = 1
scenario = yy; %edit for desired scenario

%% Parameter estimation routine lsqcurvefit
if scenario == 1 %Vm= Vma*GLU + Vmb
    coefficientsguess = [NGvalues(1:9), 1]; % initial guess for all coefficients
    LB = zeros(1,10);
elseif scenario == 2 %c_ace=c_aceA*GLU + c_aceB
    coefficientsguess = [NGvalues(1:9), 1];
    LB = zeros(1,10);
elseif scenario == 3 %c_at1=c_at1A*GLU + c_at1B
    coefficientsguess = [NGvalues(1:9), 1];
    LB = zeros(1,10);
elseif scenario == 4 %c_nonace=c_nonaceA*GLU + c_nonaceB
    coefficientsguess = [NGvalues(1:9), 1];
    LB = zeros(1,10);
elseif scenario == 5
     coefficientsguess = [NGvalues(1:9),1,1,1,1];
     LB = zeros(1,13);
else 
    coefficientsguess = NGvalues(1:9);
    LB = zeros(1,9);
end
%% Parameter estimation routine lsqcurvefit
% glucoseRASss(coefficientsguess,Xdata(1))
options = optimoptions('lsqcurvefit','TolFun',1e-16,'TolX',1e-16,'MaxFunEvals',10000);

UB = [];
if scenario == 6
    coefficients = NGHG(NGvalues);
else
    coefficients = lsqcurvefit(@(coefficients,Xdata) glucoseRASss(coefficients,Xdata,baseline,scenario,kAGT,printoutput), coefficientsguess, Xdata, Ydata,[],[],options);%,LB,UB);%,LB,UB,OPTIONS); %the function returns ANGII concentration at steady state
 

end

X = Xdata(1):0.5:Xdata(end);

for i = 1:length(X)
    Y(i)= glucoseRASss(coefficients,X(i),baseline,scenario,kAGT,printoutput);
end
for i = 1:length(Xdata)-1
    Ycalc(i)= glucoseRASss(coefficients,Xdata(i),baseline,scenario,kAGT,printoutput);
end
Ycalc(i+1)= glucoseRASss(coefficients,Xdata(end),baseline,scenario,kAGT,1);

[x,resnorm,residual,exitflag,output,lambda,jacobian] =...
         lsqcurvefit(@(coefficients,Xdata) glucoseRASss(coefficients,Xdata,baseline,scenario,kAGT,printoutput), coefficientsguess, Xdata, Ydata,[],[],options);
conf = nlparci(x,residual,'jacobian',jacobian)
%% Plots
figure(1)
hold on
plot(X,Y,'o-')
xlabel('Glucose (mM)')
% get(gca);set(gca,'FontSize',10,'FontName','Arial');
xlabel('X','FontSize',10)
ylabel('Y','FontSize',10)
% set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 3.5 2.5]);
% 
% export_fig figure(1) -r1000 -a4 -q101 -painters -eps -png -tiff

ylabel({'ANGII Change';'from Baseline'})
hold off
RMSE(yy+1) = sqrt(mean((Ydata - Ycalc).^2));
if scenario == 1 
figure(2)
subplot(221)
plot(X, coefficients(10)*X+coefficients(9))
xlabel('Glucose (mM)')
ylabel('V_m (mM/s)')
titletext = 'Scenario %i: Vm = %1.2e [Glucose] + %1.2e\n';
g= sprintf(titletext,scenario, coefficients(10),coefficients(9));
title(g)
elseif scenario == 2
figure(2)
subplot(222)
plot(X, coefficients(10)*X+coefficients(2))
xlabel('Glucose (mM)')
ylabel('c_{ace} (mM/s)')
titletext = 'Scenario %i: c_{ace} = %1.2e [Glucose] + %1.2e\n';
g= sprintf(titletext,scenario, coefficients(10),coefficients(2));
title(g)
elseif scenario == 3
figure(2)
subplot(223)
plot(X, coefficients(10)*X+coefficients(7))
xlabel('Glucose (mM)')
ylabel('c_{at1} (mM/s)')
titletext = 'Scenario %i: c_{at1} = %1.2e [Glucose] + %1.2e\n';
g= sprintf(titletext,scenario, coefficients(10),coefficients(7));
title(g)
% title(['Scenario 3: c_{at1} = ' num2str(coefficients(7)) 'Glucose + ' num2str(coefficients(9))])
elseif scenario == 4
figure(2)
subplot(224)
plot(X, coefficients(10)*X+coefficients(3))
xlabel('Glucose (mM)')
ylabel('c_{nonace} (nM/s)')
titletext = 'Scenario %i: c_{nonace} = %1.2e [Glucose] + %1.2e\n';
g= sprintf(titletext,scenario, coefficients(10),coefficients(3));
title(g)
elseif scenario == 5
figure(3)
subplot(221)
plot(X, coefficients(10)*ones(size(X)))
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

figure(1)
axis([5 40 0.8 2.2])

legend('Data from [1]',['Case 0: RMSE = ' num2str(round(RMSE(0+1),2))],...
    ['Case 1: RMSE = ' num2str(round(RMSE(1+1),2))],...
    ['Case 2: RMSE = ' num2str(round(RMSE(2+1),2))],...
    ['Case 3: RMSE = ' num2str(round(RMSE(3+1),3))],...
    ['Case 4: RMSE = ' num2str(round(RMSE(4+1),2))],...
    ['Case 5: RMSE = ' num2str(round(RMSE(5+1),8))],...
    ['Case 6: RMSE = ' num2str(round(RMSE(6+1),8))], 'Location','EastOutside')

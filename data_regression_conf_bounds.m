GLU = [5, 25, 40]';
ANG_II= [1, 1.4, 2.1]';
Ding_error = [ 0.7788094, 1.1386267, 1.305477, 2.423198]; % wrong values
x = GLU;
y = ANG_II;
E = Ding_error;
[fitresult] = fit(x,y,'power1')
prediction_interval = predint(fitresult,x,0.95); % 95% prediction interval
figure(1)
plot(fitresult,'k-')
hold on
errorbar(x,y,E,'kx')
plot2 = plot(x,prediction_interval,'k--')
hold off
stringx={'Ang II (nM)'};
stringy={'% Apoptotic Cells'};
xlabel(stringx,'fontsize',20),ylabel(stringy,'fontsize',20)
legend('fitted curve','data','prediction interval','Location','Northwest')

Jia_time = [0, 7, 14, 21, 28]';
Jia_proteinuria = [4.892368, 11.056751, 14.383562, 17.710371, 32.19178]';
Jia_error_upper = [ 2.672574, 3.913913, 5.97831, 8.432756, 11.353787]; % wrong values
Jia_error_lower = [ 1.5749786, 6.066248, 8.012159, 7.622092, 10.411316]; % wrong values
L = Jia_error_lower;
U = Jia_error_upper;
x = Jia_time;
y = Jia_proteinuria;
[fitresult,gof] = fit(x,y,'exp1')
prediction_interval = predint(fitresult,x,0.95); % 95% prediction interval
figure(2)
plot(fitresult,'k-')
hold on
errorbar(x,y,L,U,'kx')
plot2 = plot(x,prediction_interval,'k--')
hold off
stringx={'Time (days)'};
stringy={'Proteinuria (mg/day)'};
xlabel(stringx,'fontsize',20),ylabel(stringy,'fontsize',20)
legend('fitted curve','data','prediction interval','Location','Northwest')

%matlab links:
%http://www.mathworks.com/help/curvefit/list-of-library-models-for-curve-and-surface-fitting.html#btbcvnl
%http://www.mathworks.com/help/curvefit/predint.html
% http://www.mathworks.com/help/curvefit/evaluating-goodness-of-fit.html
% http://www.mathworks.com/help/curvefit/fit.html
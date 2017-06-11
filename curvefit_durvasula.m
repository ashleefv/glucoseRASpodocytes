xdata = ...
[5,25,40];
ydata = ...
[1,1.335,2.1];
fun = @(x,xdata)(x(1).*exp(x(2).*xdata)); %exp model- a*exp (b*x)
x0 = [.5,.4];
x = lsqcurvefit(fun,x0,xdata,ydata)

times = linspace(xdata(1),xdata(end));
plot(xdata,ydata,'ko',times,fun(x,times),'b-')
legend('Data','Fitted exponential')
title('Data and Fitted Curve')

[fitresult] = fit(xdata',ydata','exp1')
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[zeros(1,10)],...
               'Upper',[Inf*(ones(1,10))],...
               'StartPoint',[ones(1,10)]);
ft = fittype('glucoseRASss(GLU,coefficients)','options',fo);
[curve2,gof2] = fit(xdata,ydata,ft)
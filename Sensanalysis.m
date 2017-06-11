function Sensanalysis 
close all
format long e
%% Take data from plots

Xdata = [5];% , 25, 40]; % GLU %mmol/L
% Xdata_new = Xdata.*1.01
p = 1.1;
%% initial guesses
NG;
load('NGvalues','NGvalues')
printoutput = 1;


options = optimoptions('fmincon','Algorithm','sqp','TolFun',1e-16,'TolX',1e-16);%'MaxFunEvals',10000);
for  yy = 1:8    
event = yy;
    ANGII_newNG(yy) = glucoseRASssSens(NGvalues,Xdata,1,1,event,printoutput,p)
    ANGII_newHG(yy) = glucoseRASssSens(NGvalues,Xdata,1,2,event,printoutput,p)
end

ANGII_oldNG = glucoseRASssSens(NGvalues,Xdata,1,1,0,printoutput,p).*ones(size(ANGII_newNG));
ANGII_oldHG = glucoseRASssSens(NGvalues,Xdata,1,2,0,printoutput,p).*ones(size(ANGII_newHG))%%
%Fully normalized Sensitivity Calculation
SensAnalysisNG = abs((ANGII_oldNG-ANGII_newNG)./(ANGII_oldNG.*0.1));
SensAnalysisHG = abs((ANGII_oldHG-ANGII_newHG)./(ANGII_oldHG.*0.1));
x = 1:1:8;
 A = [SensAnalysisNG; SensAnalysisHG];
 % Sort the parameters starting from the most sensitive ones in descending
 % order
 [y,I] = sort(A(1,:),'descend');
 B=A(:,I);

figure(1);clf; 
%h=semilogy(x,B./max(SensAnalysisNG),'x');
bar(x,(B./max(SensAnalysisNG))');
    set(gca,'YScale','log');
% set(h,'linewidth',4);
% set(h,'Markersize',10);
set(gca,'Fontsize',10);
grid on;
ylabel('Sensitivity of ANG II','Fontsize',10)%,'FontWeight','Bold');
legend('NG Sensitivity','HG Sensitivity','Fontize',10,'Location','Northeast')

ax = gca;
ax.XTickLabel = {'c_{AItoII}','k_{AGT}','V_m','c_{AT1}','c_{APA}','c_{AT2}','c_{ACE2}','c_{NEP}','Fontsize',18,'FontWeight','Bold'};

set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 6 4]);
%export_fig('C:/Users/Minu/Desktop/minu_pilvankar_CHE5110/matlab/Figures/PaperVersion/SensitivityWithkAGT64', '-pdf', '-png', '-eps', '-tiff');

end


   

%Uses eFAST from Kirschner http://malthus.micro.med.umich.edu/lab/usadata/
%% PARAMETER INITIALIZATION
% set up max and mix matrices
%% initial guesses
NG;
load('NGvalues','NGvalues');
 Vm = NGvalues(1);
    kAGT = NGvalues(2);
    c_ACE = NGvalues(3); %c_ace +c_nonace
    c_nep = NGvalues(4);
    c_ace2 = NGvalues(5);
    c_apa = NGvalues(6);
    c_at1 = NGvalues(7);
    c_at2 = NGvalues(8);
scalevalue = 1e-3;

pmin=[Vm*(scalevalue), % VmaxoverKm 
kAGT*(scalevalue), % k_cat_Renin
c_ACE*(scalevalue), % k_feedback
c_nep*(scalevalue),
c_ace2*(scalevalue),
c_apa*(scalevalue),
c_at1*(scalevalue),% feedback_capacity
c_at2*(scalevalue), % k_cons_AngII
1]; % dummy

scalevalue = 1e+3;
pmax=[Vm*(scalevalue), % VmaxoverKm 
kAGT*(scalevalue), % k_cat_Renin
c_ACE*(scalevalue), % k_feedback
c_nep*(scalevalue),
c_ace2*(scalevalue),
c_apa*(scalevalue),
c_at1*(scalevalue),% feedback_capacity
c_at2*(scalevalue),
1];
% pmax=[VmaxoverKm*(scalevalue), % VmaxoverKm 
% k_cat_Renin*(scalevalue), % k_cat_Renin
% k_feedback*(scalevalue), % k_feedback
% feedback_capacity*(scalevalue), % feedback_capacity
% k_cons_AngII*(scalevalue), % k_cons_AngII
% 1]; % dummy
    % Parameter Labels 
efast_var={'k_{AGT}','c_{AItoII}','c_{RENIN}','c_{AT1}','c_{APA}','c_{AT2}','c_{ACE2}','c_{NEP}','dummy'}; %string of the labels

%             Si(i,t,u) = mean(Vi)/mean(V);
%             Sti(i,t,u) = 1-mean(Vci)/mean(V);
%             rangeSi(i,t,:,u) = Vi./V;
%             rangeSti(i,t,:,u) = 1-(Vci./V);
load('baseline','baseline');
Model_efast
% x=1:6; %5 input variables + dummy
x = 1:1:9;
% for i = x
%     for t = 1:length(Si(1,:,1))
%         for u = 1:length(Si(1,1,:))
%             meanSi(i,t,u) = mean(rangeSi(i,t,:,u));
%             stdSi(i,t,u) = std(rangeSi(i,t,:,u));
%             meanSti(i,t,u) = mean(rangeSti(i,t,:,u));
%             stdSti(i,t,u) = std(rangeSti(i,t,:,u));
%         end
%     end
% end
    
figure(1)
%b=bar(x,[Si(1:3,1,1) Sti(1:3,1,1); Si(5,1,1) Sti(5,1,1); Si(4,1,1) Sti(4,1,1); Si(6,1,1) Sti(6,1,1)]);
% [y,I] = sort(A(1,:),'descend');
%  B=A(:,I);
A = [Si(1:9,1,1) Sti(1:9,1,1)];
[y,I] = sort(A(1,:),'descend')
 B=sort(A(:,I),'descend');
 
b=bar(x,B)%[Si(1:9,1,1) Sti(1:9,1,1)]);% Si(5,1,1) Sti(5,1,1); Si(4,1,1) Sti(4,1,1); Si(6,1,1) Sti(6,1,1)]);
b(1).FaceColor = [19 106 177]/255; %blue
b(2).FaceColor = [126 162 43]/255; %green
ax = gca;
ax.XTickLabel = {efast_var{1},efast_var{2},efast_var{3},efast_var{4},efast_var{5},efast_var{6},efast_var{7},efast_var{8},efast_var{9},'FontName','Arial','FontSize',10};
% hold on
% errorbar(x,[meanSi(1:3,1,1) ; meanSi(5,1,1) ; meanSi(4,1,1); meanSi(6,1,1) ],...
%     [stdSi(1:3,1,1) ; stdSi(5,1,1) ; stdSi(4,1,1); stdSi(6,1,1) ],'.')
% hold off
legend('S_i','S_{Ti}','Location','NorthEast')
ylabel('eFAST Sensitivity of Ang II','FontName','Arial','FontSize',10)
set(gca,'FontName','Arial','FontSize',10);
axis([0 8.5 0 0.6])
%set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 3.5 2.5]);
    set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 6 4]);
% export_fig('C:/Users/Minu/Desktop/BullMathBiol/figures/globalSA', '-pdf', '-png', '-eps', '-tiff');
   
    
    

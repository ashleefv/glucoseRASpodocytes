x = [5,40];
y1 = [3.0351e-7,3.94563e-7];
get(gca);set(gca,'FontSize',10,'FontName','Arial');
yyaxis left
plot(x,y1,'r:x','LineWidth',1.5,'MarkerSize',5);
y2 = [0.016504,0.04126];
yyaxis right
plot(x,y2,'b--o','LineWidth',1.5,'MarkerSize',5)

xlabel('Glucose (mM)','FontSize',12)
ylabel('(mol/L/s)','FontSize',12)
yyaxis left
ylabel('(1/s)','FontSize',12)
legend('V_{max,renin}','c_{at1}')

set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 3.5 2.5]);

export_fig lin_vm_cat1 -r1000 -a4 -q101 -painters -eps -png -tiff
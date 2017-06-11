
function linparams = 
VmNG = coefficients(9);
VmHG = 1.3 * VmNG;
caceNG = 0.0044387;
caceHG = 1.5*caceNG;
cnonaceNG = 0.0013316;
cnonaceHG = 5*cnonaceNG;
cat1NG = 0.016504;
cat1HG = 2.5*cat1NG;
A = [5,1,0,0,0,0,0,0;...
     40,1,0,0,0,0,0,0;...
     0,0,5,1,0,0,0,0;...
     0,0,40,1,0,0,0,0;...
     0,0,0,0,5,1,0,0;...
     0,0,0,0,40,1,0,0;...
     0,0,0,0,0,0,5,1;
     0,0,0,0,0,0,40,1];
 B = [VmNG;VmHG;caceNG;caceHG;cnonaceNG;cnonaceHG;cat1NG;cat1HG];
 x = A\B
 %x = [Vma,Vmb,c_aceA,c_aceB,c_nonaceA,c_nonaceB,c_at1A,c_at1B]
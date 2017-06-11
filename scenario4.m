function output = scenario5(coefficients,GLU)
    kAGT = coefficients(1);
    c_ACEB = coefficients(2);
   % c_nonaceB = coefficients(3);
    c_nep = coefficients(3);
    c_ace2 = coefficients(4);
    c_apa = coefficients(5);
    c_at1B = coefficients(6);
    c_at2 = coefficients(7);
    Vmb= coefficients(8);
    c_ACEA = coefficients(9);
   % c_nonaceA = coefficients(11);
    c_at1A = coefficients(10);
    Vma = coefficients(11);
    Vm = Vma*GLU+Vmb;
    c_ACE = c_ACEA*GLU + c_ACEB;
   % c_nonace = c_nonaceA*GLU + c_nonaceB;
    c_at1 = c_at1A*GLU + c_at1B;
    
    output = [Vm,kAGT,c_ACE,c_nep,c_ace2,c_apa,c_at1,c_at2];
end
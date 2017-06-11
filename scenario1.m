function output = scenario1(coefficients,GLU)
    kAGT = coefficients(1);
    c_ACE = coefficients(2);
   % c_nonace = coefficients(3);
    c_nep = coefficients(3);
    c_ace2 = coefficients(4);
    c_apa = coefficients(5);
    c_at1 = coefficients(6);
    c_at2 = coefficients(7);
    Vmb= coefficients(8);
    Vma = coefficients(9);
    Vm = Vma*GLU+Vmb;
    output = [Vm,kAGT,c_ACE,c_nep,c_ace2,c_apa,c_at1,c_at2];
end
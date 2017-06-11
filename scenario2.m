function output = scenario2(coefficients,GLU)
    kAGT = coefficients(1);
    c_ACEB= coefficients(2);
    c_ACEA= coefficients(9);
    c_ACE = c_ACEA*GLU + c_ACEB;
   % c_nonace = coefficients(3);
    c_nep = coefficients(3);
    c_ace2 = coefficients(4);
    c_apa = coefficients(5);
    c_at1 = coefficients(6);
    c_at2 = coefficients(7);
    Vm= coefficients(8);
    output = [Vm,kAGT,c_ACE,c_nep,c_ace2,c_apa,c_at1,c_at2];
end

function output = scenario3(coefficients,GLU)
    kAGT = coefficients(1);
    c_ACE = coefficients(2);
    %c_nonace = coefficients(3);
    c_nep = coefficients(3);
    c_ace2 = coefficients(4);
    c_apa = coefficients(5);
    c_at1B = coefficients(6);
    c_at1A = coefficients(9);
    c_at1 = c_at1A*GLU + c_at1B;
    c_at2 = coefficients(7);
    Vm= coefficients(8);
    output = [Vm,kAGT,c_ACE,c_nep,c_ace2,c_apa,c_at1,c_at2];
end
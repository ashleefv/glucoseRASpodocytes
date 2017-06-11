function output = scenarioM1(coefficients,GLU)
    c_ace = coefficients(1);
    c_nonace = coefficients(2);
    c_at1 = coefficients(3);
    Vmb= coefficients(4);
    Vma = coefficients(5);
    Vm = Vma*GLU+Vmb;
    output = [c_ace,c_nonace,c_at1,Vm];
end


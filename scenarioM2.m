function output = scenarioM2(coefficients,GLU)
    c_aceB= coefficients(1);
    c_aceA= coefficients(2);
    c_ace = (c_aceA*GLU )^ c_aceB;
    c_nonace = coefficients(3);
    c_at1 = coefficients(4);
    Vm= coefficients(5);
    output = [c_ace,c_nonace,c_at1,Vm];
end


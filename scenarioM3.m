function output = scenarioM3(coefficients,GLU)
    
    c_ace = coefficients(1);
    c_nonace = coefficients(2);
   
    c_at1B = coefficients(3);
    c_at1A = coefficients(4);
    c_at1 = c_at1A*GLU + c_at1B;
    
    Vm= coefficients(5);
    output = [c_ace,c_nonace,c_at1,Vm];
end

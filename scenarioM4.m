function output = scenarioM4(coefficients,GLU)
   
    c_ace = coefficients(1);
    c_nonaceB = coefficients(2);
    c_nonaceA = coefficients(3);
    c_nonace = c_nonaceA*GLU + c_nonaceB;
   
    c_at1 = coefficients(4);
    
    Vm= coefficients(5);
    output = [c_ace,c_nonace,c_at1,Vm];
end
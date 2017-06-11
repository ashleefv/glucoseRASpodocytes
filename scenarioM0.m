function output = scenarioM0(coefficients)
    kAGT = coefficients(2);
    c_ACE = coefficients(3);
    c_at1 = coefficients(7);
    Vm= coefficients(1);
    
    output = [kAGT,c_ACE,c_at1,Vm];
end
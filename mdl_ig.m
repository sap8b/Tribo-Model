function y = mdl_ig(x)
    b1 = 0.019564;
    b2 = -0.027524;
    b3 = 7.5517e+07;
    
    num = b1 + b2.*x;
    den = 1.0 + b3.*x;
    
    y = num./den;
end
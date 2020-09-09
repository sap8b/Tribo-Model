function y = mdl_a_plus(x)
    b1 = 0.27422;
    b2 = 1.957;
    b3 = 4.5029;
    
    num = b1 + b2.*x;
    den = 1.0 + b3.*x;
    
    y = num./den;
end
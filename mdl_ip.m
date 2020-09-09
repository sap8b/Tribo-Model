function y = mdl_ip(x)
    b1 = 4.591e+11;
    b2 = -7.7924e+11;
    b3 = 4.3974e+14;
    
    num = b1 + b2.*x;
    den = 1.0 + b3.*x;
    
    y = num./den;
end
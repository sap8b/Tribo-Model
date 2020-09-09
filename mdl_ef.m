function y = mdl_ef(x)
    b1 = -0.29183;
    b2 = 17.015;
    b3 = 16.472;
    
    num = b1 + b2.*x;
    den = 1.0 + b3.*x;
    
    y = num./den;
end
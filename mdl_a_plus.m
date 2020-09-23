function y = mdl_a_plus(x)
    b1 = 0.38815; %0.27422;
    b2 = 2.616; %1.957;
    b3 = 11.231; %4.5029;
    
    num = b1 + b2.*x;
    den = 1.0 + b3.*x;
    
    y = num./den;
end
function y = mdl_ip(x)
    b1 = 0.026801; %4.591e+11;
    b2 = -10.946; %-7.7924e+11;
    b3 = -7.1448e-06; %4.3974e+14;
    
    y = (b1.*exp(b2.*x)) + b3;
    
%     num = b1 + b2.*x;
%     den = 1.0 + b3.*x;
%     
%     y = num./den;
end
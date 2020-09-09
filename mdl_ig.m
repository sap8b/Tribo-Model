function y = mdl_ig(x)
    b1 = 25.463; %0.019564;
    b2 = -10.695; %-0.027524;
    b3 = 0.27382; %7.5517e+07;
    
    y = (b1.*exp(b2.*x)) + b3;
%     num = b1 + b2.*x;
%     den = 1.0 + b3.*x;
%     
%     y = num./den;
end
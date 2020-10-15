function E2 = mdl_ef(U,i0_base,i0_Me)
%     b1 = -0.29183;
%     b2 = 17.015;
%     b3 = 16.472;
%     
%     num = b1 + b2.*U;
%     den = 1.0 + b3.*U;
%     
%     y = num./den;
    
    R = 8.314; %J/mol K
    T = 298; %K
    alpha = 0.5;
    z = 3;
    F = 96485; %coul/mol
    s = 0.5e-9; %m
    
    num = R*T.*U;
    den = alpha * z * F * s;
    
    curr_val = log(i0_Me./i0_base);
    
    E2 = (num./den) .* curr_val;
    
end
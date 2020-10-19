function i_me = i_me_corr(eta)
    i0 = 9.07e-16; %A/m2
    alpha_me = 0.5;
    F = 96485; %coul/mol
    R = 8.314; %J/mol K
    T = 298; %K
    
    i_me = i0 .* exp((alpha_me * F .* eta)./(R*T));
    
end
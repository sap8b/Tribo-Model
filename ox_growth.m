function node_height = ox_growth(i0g,a_p,E0,kf,eta,dt)
    F = 96485; %col/mol
    R = 8.314; %J/mol K
    T = 298; % K
    
    g_p = (a_p * F)/(R * T);
    
    first_term = 1.0/(g_p * E0);
    exp_term = exp(g_p * eta);
    inside_log = 1.0 + (g_p * E0 * kf * i0g * exp_term * dt);
    
    node_height = first_term * log(inside_log);
end
function [contact_radius_M, depth_M, p_0] = Hertzian_Contact(E_t, nu_t, E_o, nu_o, L, R)
    val_tip = (1-(nu_t^2))/E_t;
    val_oxide = (1-(nu_o^2))/E_o;
    E_star = 1.0/(val_tip + val_oxide);
    
    d = ((0.75*L)/(E_star * sqrt(R)))^(2/3);
    pi = 3.14159265;
    
    a = sqrt(R*d);
    
    p_0 = (3*L)/(2*pi*a^2);
    contact_radius_M = a;
    depth_M = d;

%     p_0 = (3.0/2.0)*L; % p_0 = p_max in GPa

%     contact_radius_M = ((pi*R)/(2*E_star)).*p_0; %m
%     depth_M = (pi/(2*E_star)).*p_0.*contact_radius_M; %m    
end
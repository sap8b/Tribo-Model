function [i,numerator,denominator] = i_growth(k_film, E_film, i0, g_plus, eta, dt)
    a_dt_modifier = 1.0; %1.3333;
    modified_dt = a_dt_modifier * dt;
    
    numerator = i0 * exp(g_plus*eta);
    denominator = 1 + (E_film * k_film * g_plus * numerator * modified_dt); %2.0* 1.87 1.0256*
    
    i = numerator/denominator;
    
    if abs(imag(numerator)) > 0.0
        disp(num)
    end    
end
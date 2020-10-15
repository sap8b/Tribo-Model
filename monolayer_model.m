function ic = monolayer_model(eta,i0_monolayer)
%     if i0_monolayer > 1.0e-30
%         k = 2.0e3; %1.0e1; %5.0e3;
%         ic = i0_monolayer.*exp(-k.*eta);        
%     else
%         ic = 0.0;
%     end
%     
    F = 96485; %coul/mol
    R = 8.314; %J/mol K
    T = 298; % K
    alpha_ox = 0.5;
    
    ic = i0_monolayer * exp((alpha_ox * F * eta)/(R*T));
    
%     ic = 0.0;

end
function im = monolayer_model(dt,i0_monolayer)
    if i0_monolayer > 1.0e-30
        k = 2.0e3; %1.0e1; %5.0e3;
        im = i0_monolayer.*exp(-k.*dt);        
    else
        im = 0.0;
    end

end
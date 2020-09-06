function ip = passive_model(dt,i0g, i0p)
    ip = zeros(1,length(dt));
    if i0g > 1.0e-30
        i0_pass = 0.0; %i0p;
    else
        i0_pass = 0.0;
    end    
    
    ip(1,:) = i0_pass;
end
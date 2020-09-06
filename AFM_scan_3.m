function [tpos, surface_heights, tot_current_density, sim_time, a_contact] = ...
    AFM_scan_3(x_pos, y_pos, z_pos, N_nodes, scan_ts, relax_time_ts, actual_ts, ...
    ts_mult, tip, sub, L, v_tip, eta_base, alpha0, i0_growth_base, E_film, i0_monolayer_base, i0_passive_base)
    %=====================================================================
    % This function models a simple AFM tip scan
    %=====================================================================
    
    %=====================================================================
    % Some simulation parameters
    %=====================================================================
    pi = 3.14159265;
%     dt = total_time/N_int;
    velocity_tip = v_tip;
    dy = y_pos(end) - y_pos(end-1);
    dx = x_pos(end) - x_pos(end-1);
    dt = ts_mult; %(dx/velocity_tip)/ts_mult;
    t_for_tip_movement = dx/velocity_tip;
    relax_ts = relax_time_ts; %round(relax_time_ts/dt);
    convert_nm_to_um = 1.0e-3;
    
    mz = min(z_pos);
    nz = max(z_pos);
    % Set the minimum oxide thickness as 10% of the minimum oxide thickness
    % in the z_pos array
    zmin1 = 0.1* min(mz);
    zmin2 = min(mz);
    zmax = max(nz);

    switch tip
        case 'Diamond'
            nu_tip = 0.10;
            E_tip = 1053e9; %Pa
            R_AFM_tip = 35.0e-9; %70.0e-9; %m
        otherwise
            nu_tip = 0.10;
            E_tip = 1053e9; %Pa    
            R_AFM_tip = 35.0e-9; %70.0e-9; %m
    end
    switch sub
        case 'Cr2O3'
            nu_substrate = 0.25; %Cr2O3 
            E_substrate = 125e9; %Pa Cr2O3 
%             H_substrate = (0.009807 * 8.25) * 1.0e9; %Pa - HV from Wikipedia for C2O3 and formula from gordonengland.co.uk
            H_substrate = (0.009807 * 50) * 1.0e9; %Pa - HV from Wikipedia for C2O3 and formula from gordonengland.co.uk            
            K_archard = 1.0e-6; %1.7e-5; %1.7e-5; % For a ferritic stainless steel in Wikipedia  
%             alpha0 = 0.32; % 0.38; % 0.41; % 0.43;
%             i0_growth_base = 1.0e-10; %10.0e-10; 50.0e-10; 100.0e-10; A/cm2
            k_film = 504; %nm.cm2/A.s 
            z = 1;
%             E_film = 0.5; %0.625; %0.581; %0.543; % 
            alpha_stress_modifier = 1;  
            rho = 5.22 *(1/(0.01 * 0.01 * 0.01)); %g/cm3 -> g/m3
            M = 151.9904; %g/mol
            Vm = M/rho; % m3/mol
        case 'UNS S32750'
            nu_substrate = 0.27; % UNS S32750 0.25; %Cr2O3 
            E_substrate = 210e9; % UNS S32750 Ps 125e9; %Pa Cr2O3 
             H_substrate = (0.009807 * 250) * 1.0e9; %Pa - HV from Wikipedia for 2507 and formula from gordonengland.co.uk
            K_archard = 5.0e-5; %1.7e-5; %1.7e-5; % For a ferritic stainless steel in Wikipedia            
              
            alpha0 = 0.25; %0.43; %0.16; %0.25; %0.25; %0.68; %0.25; %0.5;
            i0_growth_base = 1.0e-8; %A/cm2
            k_film = 504;%nm.cm2/A.s % 
            z = 1;
            E_film = 0.5;  
            alpha_stress_modifier = 1;
        otherwise
            nu_substrate = 0.2; % UNS S32750 0.25; %Cr2O3 
            E_substrate = 100e9; % UNS S32750 Ps 125e9; %Pa Cr2O3 
            H_substrate = (0.009807 * 10) * 1.0e9; %Pa - HV from Wikipedia for 2507 and formula from gordonengland.co.uk
            K_archard = 1.0e-3; %1.7e-5; % For a ferritic stainless steel in Wikipedia 
            alpha0 = 0.25; %0.68; %0.25; %0.5;
            i0_growth_base = 1.0e-10; %1.0e-16; %A
            k_film = 504; %nm.cm2/A.s   
            z = 1;
            E_film = 0.5;
            alpha_stress_modifier = 1;
    end   
    
    %=====================================================================
    % Calculate Hertzian contact parameters
    %=====================================================================
    [r_damage_m, depth_m, p_max] = Hertzian_Contact(E_tip, nu_tip, E_substrate, nu_substrate, L, R_AFM_tip);
    r_damage_nm = r_damage_m * 1.0e9;
        
    if r_damage_nm > (R_AFM_tip * 1.0e9)
       r_damage_nm = R_AFM_tip * 1.0e9;
    end
    %=====================================================================
    
    %=====================================================================
    % alpha model
    %=====================================================================
    p_app_max = 1.30409e10; %6.0529e10; %2.4688e-4; % nm for L_max = 10.0 N (K_archard/H_substrate)*p_max*dt*(velocity_tip*1.0e-9); %m 
    p_app_min = 0.0; %1.0e10;
    alpha_min = alpha0; %0.2*alpha_max_current; %0.9767*alpha_max_current; %0.42;
    alpha_max = 1.01 * alpha0; %0.43;
    p_ratio = p_max/p_app_max;
    alpha_max_current = alpha_max * (1 + p_ratio);
    
    a_slope = (alpha_max - alpha_min)/(p_app_max - p_app_min);
    %=====================================================================
    %=====================================================================
    % i0 growth model - initial damage
%     i0_slope = (i0_max - i0_min)/(p_app_max - p_app_min);
    i0_model = i0_growth_base; %i0_pre * exp(i0_exp * L);
    
    % i0 growth model - post damage
    %=====================================================================
    
%     check_y_nodes_shift = round((R_AFM_tip * 1.0e9)/dy);
    check_y_nodes_shift = 1; 
    
    if check_y_nodes_shift < 1
        num_y_nodes_shift = 1;
    else
%         num_y_nodes_shift = 1; 
        num_y_nodes_shift = check_y_nodes_shift;
    end
    
    depth_nm = depth_m * 1.0e9;
    
    Faraday_Constant = 96485; %coul/mol
    R = 8.314; %J/mol K
    T = 298; %K        
    
    %=====================================================================
    %=====================================================================
    % Create an oxide sructure to keep track of damaged and undamaged areas
    %=====================================================================
    node_counter = 1;
    for i = 1:(N_nodes(2))
        for j = 1:(N_nodes(1))
            oxide(node_counter).grid = [x_pos(j), y_pos(i)];
            oxide(node_counter).nodes = [j,i];
            oxide(node_counter).height = z_pos(j,i);
            oxide(node_counter).base_height = z_pos(j,i); %zmax; %2.5; %
            oxide(node_counter).i0_growth = 0.0;
            oxide(node_counter).alpha_node = 0.0;
            oxide(node_counter).k_film = k_film;
            oxide(node_counter).num = 0.0;
            oxide(node_counter).denom = 0.0;
            oxide(node_counter).damage_current = 0.0;
            oxide(node_counter).initiation_time = 0.0;
            oxide(node_counter).rebuild_time = 0.0;
            oxide(node_counter).has_damage = 0;
            node_counter = node_counter + 1;
        end
    end
    num_nodes = node_counter - 1;
    surface_heights = zeros(num_nodes,1);
    %=====================================================================
    
    %=====================================================================
    % Start the time evolution
    %=====================================================================    
    tpos = zeros(scan_ts+relax_ts+1,4);
    tot_current_density = zeros(scan_ts+relax_ts+1,1);
    sim_time = zeros(scan_ts+relax_ts+1,1);
    
    tip_pos_x = 0.0; %0.0:(v_tip*dt):(v_tip * total_time);
    tip_pos_y = 0.0; %(ymax/2.0) - dy; %tip_pos_x;
    direction_of_x_travel = 1; % 1 = forward, -1 = backward    
    direction_of_y_travel = 1; % 1 = up, -1 = bacvkward
    total_pause = 1; %10; 
    pause_counter = 1;
    %=====================================================================
    %=====================================================================
    % Start the area calculation
    area_sum = 0.0;
    %=====================================================================
    %=====================================================================
    % Time iterations for AFM scanning and relaxatgion
    %=====================================================================
    tip_ts = 1;
    for idx_time = 1:actual_ts%(scan_ts+relax_ts+1)
       sim_time(idx_time,1) = (idx_time-1)*dt;
       
%        if mod(idx_time,ts_mult) == 0
%            tip_ts = tip_ts + 1;
%        end
       if tip_ts >= (scan_ts+1)
           velocity_tip = 0.0;
       end
       if tip_ts >= (scan_ts+relax_ts+1)
           break;
       end
       
       % Check to see if the y-position of the AFM tip has reached the
       % upper end of the computational cell, and, if so, stop it
       if direction_of_y_travel > 0
           if tip_pos_y >= y_pos(end)
               if pause_counter <= total_pause
                   direction_of_x_travel = 0;
%                    direction_of_y_travel = 0;
                   velocity_tip = 0.0;
                   pause_counter = pause_counter + 1;
               else
                   direction_of_x_travel = 1;
                   direction_of_y_travel = -1;
                   velocity_tip = v_tip;
                   pause_counter = 1;                
               end               
           end           
       elseif direction_of_y_travel < 0
           if tip_pos_y <= 0.0
               if pause_counter <= total_pause
                   direction_of_x_travel = 0;
%                    direction_of_y_travel = 0;
                   velocity_tip = 0.0;
                   pause_counter = pause_counter + 1;
               else
                   direction_of_x_travel = 1;
                   direction_of_y_travel = 1;
                   velocity_tip = v_tip;
                   pause_counter = 1;  
               end               
           end             
       end
       
       % Determine the x,y-position of the AFM tip and if it is scanning
       % forward or backward
       if direction_of_y_travel > 0
           if tip_pos_x >= x_pos(end) && direction_of_x_travel > 0
               direction_of_x_travel = -1;
               tip_pos_y = tip_pos_y + (num_y_nodes_shift * dy);
           elseif tip_pos_x <= (0.0) && direction_of_x_travel < 0
               direction_of_x_travel = 1;
               tip_pos_y = tip_pos_y + (num_y_nodes_shift * dy);
           end           
       elseif direction_of_y_travel < 0
           if tip_pos_x >= x_pos(end) && direction_of_x_travel > 0
               direction_of_x_travel = -1;
               tip_pos_y = tip_pos_y - (num_y_nodes_shift * dy);
           elseif tip_pos_x <= (0.0) && direction_of_x_travel < 0
               direction_of_x_travel = 1;
               tip_pos_y = tip_pos_y - (num_y_nodes_shift * dy);
           end           
       end
       
       if velocity_tip > 0.0
           %=====================================================================
           % AFM tip scanning
           %=====================================================================
           for j = 1:num_nodes 
                testy = oxide(j).grid(2);
                testx = oxide(j).grid(1);
                testr = sqrt((testx - tip_pos_x)^2 + (testy - tip_pos_y)^2);
                check_node_height = oxide(j).height;
                check_damage = oxide(j).has_damage;

                node_base_height = oxide(j).base_height;
                
                if (testr < r_damage_nm) 
                    % Checks to see if the test node is inside the damage
                    % area
                    if  (check_damage < 0.5) && ( check_node_height > zmin1)
                        % The test node is undamaged but is inside the
                        % damage area
                        oxide(j).has_damage = 1;
                        init_time = sim_time(idx_time,1);
                        oxide(j).initiation_time = init_time;
                        oxide(j).rebuild_time = init_time; % + dt

                        delta_t = sim_time(idx_time,1) - oxide(j).initiation_time;
                        
                        ar_ratio = (testr^2)/(r_damage_nm^2);
                        p_node = p_max * sqrt(1 - (ar_ratio));
                        delta_h_m = (K_archard/H_substrate)*p_node*dt*(velocity_tip*1.0e-9); %m
                        delta_h_nm = delta_h_m * 1.0e9; %nm
                        oxide(j).height = check_node_height - delta_h_nm;                
                        node_height_nm = oxide(j).height; % nm                        
                        
                        height_ratio = node_height_nm/node_base_height;
                        eta_adj = (p_node/E_substrate)*(2*p_node*Vm)/(3*3*Faraday_Constant);
                        if height_ratio < 1.0
                            
                            
                            alpha = alpha0; %(a_slope * p_node) + alpha_min; %alpha0 * eta_adj; %*(1+alpha_stress_modifier*(p_node/E_substrate));
                            if alpha > alpha_max
                                disp(alpha)
                            end                        
                            i0_growth = i0_model; %i0_slope * p_node + i0_min;  
                            i0_monolayer = i0_monolayer_base;
                        else
                            alpha = alpha0;
                            eta_adj = 0.0;
                            i0_growth = 0.0;
                            i0_monolayer = 0.0;
                        end
                            
                    elseif  (check_damage > 0.5) && ( check_node_height > zmin1)

                        % The test node has already been damaged and is
                        % inside the damage area
                        
                        init_time = oxide(j).initiation_time;
                        delta_t = sim_time(idx_time,1) - init_time;
                        
                        ar_ratio = (testr^2)/(r_damage_nm^2);
                        p_node = p_max * sqrt(1 - (ar_ratio));
                        delta_h_m = (K_archard/H_substrate)*p_node*dt*(velocity_tip*1.0e-9); %m
                        delta_h_nm = delta_h_m * 1.0e9; %nm
                        oxide(j).height = check_node_height - delta_h_nm;                
                        node_height_nm = oxide(j).height; % nm                        

                        height_ratio = node_height_nm/node_base_height;
                        eta_adj = (p_node/E_substrate)*(2*p_node*Vm)/(3*3*Faraday_Constant);
                        
                        if height_ratio < 1.0
                            
                            alpha = alpha0; %(a_slope * p_node) + alpha_min; %
                            if alpha > alpha_max
                                disp(alpha)
                            end                           
                            i0_growth = i0_model; %i0_slope * p_node + i0_min;  
                            i0_monolayer = i0_monolayer_base;
                        else
                            alpha = alpha0;
                            eta_adj = 0.0;
                            i0_growth = 0.0;
                            i0_monolayer = 0.0;
                        end
%                         continue

                    elseif  (check_damage < 0.5) && ( check_node_height <= zmin1)
                        % Test node is undamaged but the asperity height is
                        % zero or negative so just zero it out
                        oxide(j).has_damage = 0.0;
                        oxide(j).height = zmin1; 
                        surface_heights(j) = oxide(j).height;
                        continue

                    elseif  (check_damage > 0.5) && ( check_node_height <= zmin1)

                        % Test node has damage that has reduced its height
                        % to zmin or below. Set height to zmin and keep
                        % monitoring the corrosion current
                        oxide(j).height = zmin1; %0.0;
                        node_height_nm = oxide(j).height;
                        
                        init_time = oxide(j).initiation_time;
                        delta_t = sim_time(idx_time,1) - init_time;

                        height_ratio = node_height_nm/node_base_height;
                        eta_adj = (p_node/E_substrate)*(2*p_node*Vm)/(3*3*Faraday_Constant);
                        
%                         delta_h_nm = oxide(j).base_height - node_height_nm;
                        if height_ratio < 1.0
                           
                            alpha = alpha0; %(a_slope * p_node) + alpha_min; %
                            if alpha > alpha_max
                                disp(alpha)
                            end                            
                            
                            i0_growth = i0_model; %i0_growth_base; %i0_slope * p_node + i0_min;
                            i0_monolayer = i0_monolayer_base;
                        else
                            alpha = alpha0;
                            eta_adj = 0.0;
                            i0_growth = 0.0;
                            i0_monolayer = 0.0;
                        end
%                         continue

                    end

                else
                    % The test node is outside the damage area from the AFM
                    % tip but, if it was previously damaged, keep
                    % monitoring the oxide re-formation current
                    if check_damage > 0.5
                        init_time = oxide(j).initiation_time;
                        delta_t = sim_time(idx_time,1) - init_time;
                        
                        node_height_nm = oxide(j).height;
                        height_ratio = node_height_nm/node_base_height;
%                         delta_h_nm = oxide(j).base_height - node_height_nm;                        
                        if height_ratio < 1.0
                            eta_adj = 0.0; %1.0 - 
                            
                            alpha = alpha0; %alpha_min; %alpha_max_current; %alpha0; %(a_slope * delta_h_nm) + alpha_min; %
                            if alpha > alpha_max
                                disp(alpha)
                            end                            
                            i0_growth = i0_model; %i0_slope * p_node + i0_min;  
                            i0_monolayer = i0_monolayer_base;
                        else
                            alpha = alpha0;
                            eta_adj = 0.0;
                            i0_growth = 0.0;
                            i0_monolayer = 0.0;
                        end
                        
                    else
                        % Node is outside the damage area and is undamaged
                        i0_growth = 0.0;
                        alpha = alpha0;
                        i0_monolayer = 0.0;
                    end
                    
                end                 
                
                g_plus = (alpha * z * Faraday_Constant)/ (R * T); 
                %=====================================================================
                eta_adj = 0.0;
                %=====================================================================
                eta = eta_base + eta_adj;
                [temp__interface_current,num,denom] = i_growth(k_film, E_film, i0_growth, g_plus, eta, delta_t);
                temp_mono_current = monolayer_model(delta_t,i0_monolayer);
                temp_passs_current = passive_model(delta_t, i0_growth, i0_passive_base);

                oxide(j).damage_current = temp__interface_current + temp_mono_current + temp_passs_current; 
                oxide(j).i0_growth = i0_growth;
                oxide(j).alpha_node = alpha;
                oxide(j).num = num;
                oxide(j).denom = denom;                

                surface_heights(j) = oxide(j).height;
               
           end

       else
           %=====================================================================
           % AFM tip stopped - monitor relaxation current
           %=====================================================================
           for j = 1:num_nodes
               if oxide(j).has_damage > 0.5
                    init_time = oxide(j).initiation_time;
                    delta_t = sim_time(idx_time,1) - init_time;
                    
                    node_height_nm = oxide(j).height;% * 1.0e-3; nm                    
                    height_ratio = node_height_nm/oxide(j).base_height;
%                     delta_h_nm = oxide(j).base_height - node_height_nm; 
                    if height_ratio < 1.0
                        eta_adj = 0.0; %1.0 - height_ratio;

                        alpha = alpha0; %alpha_min; %alpha_max_current; %alpha0; %(a_slope * delta_h_nm) + alpha_min; %
                        if alpha > alpha_max
                            disp(alpha)
                        end                        
                        
                        i0_growth = i0_growth_base; %i0_slope_2 * height_ratio + i0_max;
                        i0_monolayer = i0_monolayer_base;
                    else
                        alpha = alpha0;
                        
                        eta_adj = 0.0;
                        i0_growth = 0.0;
                        i0_monolayer = 0.0;
                    end
                    
                    g_plus = (alpha * z * Faraday_Constant)/ (R * T);                   
                    [temp__interface_current,num,denom] = i_growth(k_film, E_film, i0_growth, g_plus, eta, delta_t);
                    temp_mono_current = monolayer_model(delta_t,i0_monolayer);
                    temp_passs_current = passive_model(delta_t, i0_growth, i0_passive_base);

                    oxide(j).damage_current = temp__interface_current + temp_mono_current + temp_passs_current; 
                    oxide(j).i0_growth = i0_growth;
                    oxide(j).alpha_node = alpha;
                    oxide(j).num = num;
                    oxide(j).denom = denom;                

                    surface_heights(j) = oxide(j).height;                  
               end             
               surface_heights(j) = oxide(j).height;
           end     
       end           
        
       tpos(idx_time,1) = tip_pos_x;
       tpos(idx_time,2) = tip_pos_y;
       tpos(idx_time,3) = tip_pos_y + r_damage_nm;
       tpos(idx_time,4) = tip_pos_y - r_damage_nm;
       
       tip_pos_x = tip_pos_x + (direction_of_x_travel * dx);
       area_sum = area_sum + (((2*r_damage_nm) * dy)*1.0e-7*1.0e-7);
       tip_ts = tip_ts + 1;
%        if mod(idx_time,ts_mult) == 0
% 
%        end              
       
       %=====================================================================
       % Sum total damage current
       %=====================================================================       
       temp_tot_current = 0.0;
       for j = 1:num_nodes
           i_t = oxide(j).damage_current; %* dx*1.0e-7 * dy*1.0e-7
           temp_tot_current = temp_tot_current + i_t; %converted to A           
       end
       tot_current_density(idx_time,1) = temp_tot_current/num_nodes;
       a_contact = area_sum;
       
       %=====================================================================
       % Rebuild oxide thickness
       %=====================================================================
%        for j = 1:num_nodes
%            if oxide(j).has_damage > 0.5
%                
%                delta_t = sim_time(idx_time,1) - oxide(j).rebuild_time;
%                
%                g_plus = (oxide(j).alpha_node * z * Faraday_Constant)/ (R * T);
%                oxide_dH = (1.0/(g_plus*E_film)) * log(oxide(j).denom);
%                
%                old_height = oxide(j).height;
%                new_height = old_height + oxide_dH;
%                
%                oxide(j).height = new_height; 
%                
%            end
% 
%        end
        
    end
    %=====================================================================

end
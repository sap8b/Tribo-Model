function [tpos, surface_heights, tot_current_density, sim_time, a_contact] = ...
    AFM_scan_4(x_pos, y_pos, z_pos, N_nodes, scan_ts, relax_time_ts, actual_ts, ...
    dt, tip, sub, L, v_tip, eapp_base, ecorr_base, alpha0, i0_growth_base, E2_0, ...
    i0, i0Me, i0_monolayer_base, i0_passive_base, cutoff, v_act)
    %=====================================================================
    % This function models a simple AFM tip scan
    %=====================================================================
    % Define a few physical constants
    %=====================================================================
    Faraday_Constant = 96485; %coul/mol
    R = 8.314; %J/mol K
    T = 298; %K    
    kb = 1.38e-23;
    %=====================================================================
    % Define some new variables from the parameters passed to the function
    %=====================================================================
    velocity_tip = v_tip;
    dy = y_pos(end) - y_pos(end-1);
    dx = x_pos(end) - x_pos(end-1);
%     relax_ts = relax_time_ts;     
%     mz = min(z_pos);
    % Set the minimum oxide thickness as 10% of the minimum oxide thickness
    % in the z_pos array
%     zmin1 = 0.1* min(mz);
%     b_cutoff = 0.01; %0.25;  
    another_eta_adjuster = 0.15; %1.0; %0.9;
    %=====================================================================
    %=====================================================================
    % Define material parameters for different simulation cases
    %=====================================================================
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
            H_substrate = (0.009807 * 8.25) * 1.0e9; %Pa - HV from Wikipedia for Cr2O3 and formula from gordonengland.co.uk            
            K_archard = 3.0e-5; %3.5e-5; % 5.0e-5; %5.72e-5; %1.7e-5; % For a ferritic stainless steel in Wikipedia  
            k_film = 504; %nm.cm2/A.s 
            z = 1;
            alpha_stress_modifier = 1;  
            rho = 5.22 *(1/(0.01 * 0.01 * 0.01)); %g/cm3 -> g/m3
            M = 151.9904; %g/mol
            Vm = M/rho; % m3/mol
        case 'UNS S32750'
            nu_substrate = 0.27; % UNS S32750 0.25; %Cr2O3 
            E_substrate = 210e9; % UNS S32750 Ps 125e9; %Pa Cr2O3 
             H_substrate = (0.009807 * 250) * 1.0e9; %Pa - HV from Wikipedia for 2507 and formula from gordonengland.co.uk
            K_archard = 5.0e-5; %1.7e-5; %1.7e-5; % For a ferritic stainless steel in Wikipedia            
            k_film = 504;%nm.cm2/A.s % 
            z = 1;
            alpha_stress_modifier = 1;
        otherwise
            nu_substrate = 0.2; % UNS S32750 0.25; %Cr2O3 
            E_substrate = 100e9; % UNS S32750 Ps 125e9; %Pa Cr2O3 
            H_substrate = (0.009807 * 10) * 1.0e9; %Pa - HV from Wikipedia for 2507 and formula from gordonengland.co.uk
            K_archard = 1.0e-3; %1.7e-5; % For a ferritic stainless steel in Wikipedia 
            k_film = 504; %nm.cm2/A.s   
            z = 1;
            alpha_stress_modifier = 1;
    end   
    
    %=====================================================================
    % Calculate Hertzian contact parameters
    %=====================================================================
    [r_damage_m, depth_m, p_max] = Hertzian_Contact(E_tip, nu_tip, E_substrate, nu_substrate, L, R_AFM_tip);
    r_damage_nm = r_damage_m * 1.0e9;
%     depth_nm = depth_m * 1.0e9;
        
    % This conditional check puts limits the damage radius to the AFM tip
    % radius.  This restriction may be lifted after further testing
    if r_damage_nm > (R_AFM_tip * 1.0e9)
       r_damage_nm = R_AFM_tip * 1.0e9;
    end
    %=====================================================================
    
    %=====================================================================
    % This section to be used to model changes in $\alpha^{+}$ as a
    % function of stress
    %=====================================================================

    %=====================================================================
    % This section to be used to model changes in $i_{0,field}$ as a
    % function of stress
    %=====================================================================
%     i0_model = i0_growth_base; 
    %=====================================================================
    % Determine the number of y nodes to shift during scanning. Currently,
    % the simulation provides realistic values only if the y-shift is 1
    %=====================================================================
    check_y_nodes_shift = round((R_AFM_tip * 1.0e9)/dy);    
    if check_y_nodes_shift < 1
        num_y_nodes_shift = 1;
    else
        num_y_nodes_shift = check_y_nodes_shift;
    end
    
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
            oxide(node_counter).rebuild_height = 0.0;
            oxide(node_counter).base_rebuild_height = z_pos(j,i); %(0.8691 * (p_max/1.259e10)^-0.1827)*z_pos(j,i); %
            oxide(node_counter).old_current = 0.0;
            oxide(node_counter).new_current = 0.0;
            oxide(node_counter).i0_growth = i0_growth_base;
            oxide(node_counter).alpha_node = 0.0;
            oxide(node_counter).k_film = k_film;
            oxide(node_counter).num = 0.0;
            oxide(node_counter).denom = 0.0;
            oxide(node_counter).damage_current = 0.0;
            oxide(node_counter).initiation_time = 0.0;
            oxide(node_counter).rebuild_time = 0.0;
            oxide(node_counter).has_damage = 0;
            oxide(node_counter).cutoff_state = 0; %Off = 0, On = 1
            oxide(node_counter).cutoff_current = 0.0;
            oxide(node_counter).cutoff_time = 0.0;
            node_counter = node_counter + 1;
        end
    end
    num_nodes = node_counter - 1;
    surface_heights = zeros(num_nodes,1);
    %=====================================================================
    
    %=====================================================================
    % Start time evolution
    %=====================================================================    
    tpos = zeros(actual_ts,4);
    tot_current_density = zeros(actual_ts,1);
    sim_time = zeros(actual_ts,1);
    
    tip_pos_x = 0.0; %0.0:(v_tip*dt):(v_tip * total_time);
    tip_pos_y = 0.0; %(ymax/2.0) - dy; %tip_pos_x;
    direction_of_x_travel = 1; % 1 = forward, -1 = backward    
    direction_of_y_travel = 1; % 1 = up, -1 = bacvkward
    total_pause = 1; %10; 
    pause_counter = 1;    
    %=====================================================================
    %=====================================================================
    % Initialize the damaged area calculation
    area_sum = 0.0;
    %=====================================================================
    %=====================================================================
    % Time iterations for AFM scanning and relaxation
    %=====================================================================
    tip_ts = 1;
    for idx_time = 1:actual_ts
       sim_time(idx_time,1) = (idx_time-1)*dt; %s

       if tip_ts >= (scan_ts+1)
           % Stop the AFM tip once the scanning duration reached
           velocity_tip = 0.0;
       end
       
       % This conditional only necessary if the total number of time-steps
       % allows for forward and reverse scans of the AFM tip
       if direction_of_y_travel > 0
           if tip_pos_y >= y_pos(end)
               % AFM tip has reached the upper edge of the computational
               % cell
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
               % AFM tip has reached the lower edge of the computational
               % cell               
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
                
                %=====================================================================
                % Check for 3 conditions, the 4th can't occur
                %=====================================================================
                if (testr > r_damage_nm) && (check_damage < 0.5)
                    % Node outside tip contact and has not been damaged
                    temp__interface_current = 0.0;
                    temp_pass_current = 0.0;
                    eta_me = 0.0;
                    temp_mono_current = 0.0; %monolayer_model(eta_me,i0Me);
                    oxide(j).damage_current = temp__interface_current + temp_mono_current + temp_pass_current;
                    continue;
                    
                elseif (testr < r_damage_nm) && (check_damage < 0.5)
                    % Node inside tip contact, and is damaged in this
                    % time-step
                    oxide(j).has_damage = 1;
                    init_time = sim_time(idx_time,1);
                    oxide(j).initiation_time = init_time;
                    
                    ar_ratio = (testr^2)/(r_damage_nm^2);
                    p_node = p_max * sqrt(1 - (ar_ratio));                    
                    delta_t = sim_time(idx_time,1) - oxide(j).initiation_time;

                    delta_h_m = (K_archard/H_substrate)*p_node*dt*(velocity_tip*1.0e-9); %m
                    delta_h_nm = delta_h_m * 1.0e9; %nm

                    oxide(j).height = check_node_height - delta_h_nm;                
                    node_height_nm = oxide(j).height; % nm
                    oxide(j).rebuild_height = node_height_nm;
                    E2 = abs(mdl_ef(node_height_nm*1.0e-9, i0, i0Me));
                    
                    eta_modifier = another_eta_adjuster*(p_node/1.259e10); %E_substrate 
                    eta_adj = eta_modifier *(2*p_node*Vm)/(3*z*Faraday_Constant); %

                    eta_base1 = eapp_base - ecorr_base;
                    eta_base2 = eta_base1 - E2_0;
                    eta_me = eta_adj;  %0.0; %eta_base2 - eta_adj; % eta_adj; %
                    
                    g_plus = (alpha0 * Faraday_Constant)/ (R * T); 
                    stress_effect = 0.0; %(p_node * v_act)/(kb*T);
                    eta = (E2_0 * oxide(j).base_height) - (E2 * node_height_nm);  

                    [temp__interface_current,num,denom] = i_growth(k_film, E2_0, oxide(j).i0_growth, g_plus, stress_effect, eta + eta_me, delta_t); %
                    temp_mono_current = 0.0; %monolayer_model(eta_me,i0Me);
                    temp_pass_current = 0.0;

                    oxide(j).old_current = 0.0;
                    oxide(j).new_current = temp__interface_current;  

                    oxide(j).damage_current = temp__interface_current + temp_mono_current + temp_pass_current;
                        
%                 elseif (testr < r_damage_nm) && (check_damage > 0.5)
%                     % Node inside the tip contact and damage happened in a
%                     % previous time-step
%                     disp('Can this condition happen?')
                    
                elseif (testr > r_damage_nm) && (check_damage > 0.5)
                    % Node outside the tip contact and damage happened in a
                    % previous time-step
                    init_time = oxide(j).initiation_time;
                    delta_t = sim_time(idx_time,1) - init_time;

                    node_height_nm = oxide(j).height;
%                     height_ratio = node_height_nm/node_base_height;

                    eta_modifier = another_eta_adjuster*(p_node/1.259e10); %E_substrate 
                    eta_adj = eta_modifier *(2*p_node*Vm)/(3*z*Faraday_Constant); %
% 
%                     eta_base1 = eapp_base - ecorr_base;
%                     eta_base2 = eta_base1 - E2_0;
                    eta_me = 0.0; %eta_base2 -  %                         

                    E2 = abs(mdl_ef(node_height_nm*1.0e-9, i0, i0Me));
                    
                    g_plus = (alpha0 * Faraday_Constant)/ (R * T); 
                    stress_effect = 0.0; %(p_node * v_act)/(kb*T);
                    eta = (E2_0 * oxide(j).base_height) - (E2 * node_height_nm);  

                    [temp__interface_current,num,denom] = i_growth(k_film, E2_0, oxide(j).i0_growth, g_plus, stress_effect, eta + eta_me, delta_t); %
                    temp_mono_current = 0.0; %monolayer_model(eta_me,i0Me);
                    temp_pass_current = 0.0;  

                    oxide(j).old_current = oxide(j).new_current;
                    oxide(j).new_current = temp__interface_current; 
                    
                    oxide(j).damage_current = temp__interface_current + temp_mono_current + temp_pass_current;
                end
           end

       else
           %=====================================================================
           % AFM tip stopped - monitor relaxation current
           %=====================================================================
           for j = 1:num_nodes
               % Check 2 conditions
               if oxide(j).has_damage < 0.5                   
                   % Node no longer damaged
                   temp__interface_current = 0.0;
                   temp_pass_current = 0.0;
                   eta_me = 0.0;
                   temp_mono_current = monolayer_model(eta_me,i0Me);
                   oxide(j).damage_current = temp__interface_current + temp_mono_current + temp_pass_current;
                   continue;     
                   
               elseif oxide(j).has_damage > 0.5
                   
                   % Node still daamaged, monitor current
                   init_time = oxide(j).initiation_time;
                   delta_t = sim_time(idx_time,1) - init_time;

                   node_height_nm = oxide(j).height;

                   eta_me = 0.0; %eta_base2 - eta_adj; %                         

                   E2 = abs(mdl_ef(node_height_nm*1.0e-9, i0, i0Me));
                    
                   g_plus = (alpha0 * Faraday_Constant)/ (R * T); 
                   stress_effect = 0.0; %(p_node * v_act)/(kb*T);
                   eta = (E2_0 * oxide(j).base_height) - (E2 * node_height_nm);  

                   [temp__interface_current,num,denom] = i_growth(k_film, E2_0, oxide(j).i0_growth, g_plus, stress_effect, eta + eta_me, delta_t); %
                   temp_mono_current = 0.0; %monolayer_model(eta_me,i0Me);
                   temp_pass_current = 0.0;  

                   oxide(j).old_current = oxide(j).new_current;
                   oxide(j).new_current = temp__interface_current; 
                    
                   oxide(j).damage_current = temp__interface_current + temp_mono_current + temp_pass_current; 
               end
               
           end     
       end           
        
       tpos(idx_time,1) = tip_pos_x;
       tpos(idx_time,2) = tip_pos_y;
       tpos(idx_time,3) = tip_pos_y + r_damage_nm;
       tpos(idx_time,4) = tip_pos_y - r_damage_nm;
       
       tip_pos_x = tip_pos_x + (direction_of_x_travel * dx);
       area_sum = area_sum + (((2*r_damage_nm) * dy)*1.0e-7*1.0e-7); %cm2
       tip_ts = tip_ts + 1;
       
       %=====================================================================
       % Sum current from all nodes and rebuild the damaged oxide film by
       % calculating the overpotential for oxide formation and the amount
       % of oxide added per time-step
       %=====================================================================       
       temp_tot_current = 0.0;
       
       for j = 1:num_nodes
           i_t = oxide(j).damage_current; 
           
           if oxide(j).rebuild_height < oxide(j).base_rebuild_height  && oxide(j).has_damage > 0.5
               add_height = ((oxide(j).new_current + oxide(j).old_current)/2.0) * dt * k_film;
               new_height = oxide(j).rebuild_height + add_height;
               oxide(j).rebuild_height = new_height;
               
%                if j == 1 %&& sim_time(idx_time,1) > 60
%                    d_string = strcat('t = ',num2str(sim_time(idx_time,1)), ' Base height = ', num2str(oxide(j).base_rebuild_height), ',',' Current height = ', num2str(oxide(j).rebuild_height));
%                    disp(d_string)
%                end
           elseif oxide(j).rebuild_height >= oxide(j).base_rebuild_height  && oxide(j).has_damage > 0.5
%                new_height = 0.0;
               oxide(j).height = oxide(j).base_height;
               oxide(j).has_damage = 0.0;
           end
           
           temp_tot_current = temp_tot_current + i_t; %converted to A           
       end
       tot_current_density(idx_time,1) = temp_tot_current/num_nodes;
       a_contact = area_sum;
    end
    %=====================================================================

end
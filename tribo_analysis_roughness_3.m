function tribo_analysis_roughness_3(total_exp_time, time_scanning, ...
            r_time, exp_files, output_files, ...
            L_applied, eta, surface_morphology, a0, i0, E_f, i0m, i0p, ...
            surface_dimensions, delta_pos, nodes, v_tip, fig_num)
    %==========================================================================
    % Simulation of tribocorrosion damage
    %========================================================================== 
    data = csvread(exp_files);

    x_pos = 0.0:delta_pos(1):surface_dimensions(1);
    y_pos = 0.0:delta_pos(2):surface_dimensions(2); 
    
    %=====================================================================
    % Set-up initial surface morphology and plot it
    %=====================================================================
    [X,Y] = meshgrid(y_pos, x_pos);
    num_asperities = 100;    
    
    initial_surface_roughness = build_tribo_surface(surface_morphology, num_asperities, x_pos, y_pos, surface_dimensions(3), surface_dimensions(4));
    
    maz = max(initial_surface_roughness);
    miz = min(initial_surface_roughness);
    peak_height = max(maz);
    peak_valley = min(miz);
    delta_height = (peak_height - peak_valley)/15;
    
    num_levels = peak_valley:delta_height:peak_height; %15;
    
    plotno = 1;
    tp_true = 0;
    tpos = 0.0;
%     plot_contour_tribo(plotno,X,Y,initial_surface_roughness,num_levels,tp_true,tpos)
    %=====================================================================
    %=====================================================================
    % Set up the AFM tip scanning velocity, the number of time steps to
    % take for scanning, for relazation, and dt for scanning
    %===================================================================== 
    sliding_time_per_scan = time_scanning; %[61.0, 90.0, 140.0, 142.0, 75.0, 75.0]; %60; %120.0; %s
    relaxation_time = r_time; %total_time_per_scan(6) - sliding_time_per_scan(6);
    
    number_passes = 1;
    node1 = nodes(1);
    node2 = nodes(2);
    dx = delta_pos(1);
    
    dt_scan = dx/v_tip;   
    actual_dt = dt_scan;
    
    number_of_intervals = number_passes*round(sliding_time_per_scan/actual_dt); %round(nodes(1) * nodes(2)); %(round (1.6*(nodes(1) * nodes(2))) + 10); %        
    relaxation_time_intervals = round(relaxation_time/actual_dt);
    

    number_of_time_steps = number_of_intervals+relaxation_time_intervals+1;
    
%     %=====================================================================
%     % Set up the tip-to-oxide interaction parameters
%     %=====================================================================
%     
%     L_base = [2.0, 3.0, 4.0, 5.0, 0.009, 0.03, 0.0495, 0.064, 0.09]; %2.0:1.0:10.0;
%     L_applied = L_base.*1.0e-5; % N(?) Pa
%     %=====================================================================
%     %=====================================================================
%     % Setup the electrochemical system parameters
%     %=====================================================================    
%     Eapp = [0.0, 0.2, 0.4, 0.54]; %[0.2, 0.6]; %VSCE  0.54
%     Ecorr = -0.3; %VSCE
%     eta = Eapp - Ecorr; %V

    % Perform the tribometer scan
    total_runs = 1; %(num_potentials * num_loads) + num_a + num_i0 + num_Ef;
    i0p2 = i0p/(node1 * node2);
    
    [tpos, final_surface_roughness, i_tot_all, t_all, contact_area] = ...
        AFM_scan_3(x_pos, y_pos, initial_surface_roughness, nodes, ...
        number_of_intervals, relaxation_time_intervals, number_of_time_steps, actual_dt, ...
        'Diamond', 'Cr2O3', L_applied, v_tip, eta, a0, i0, E_f, i0m, i0p2);

    %=====================================================================
    % Update the z_vals of the asperities for plotting
    %=====================================================================    
%     node_counter = 1;
%     for i = 1:(nodes(1)+1)
%         for j = 1:(nodes(2)+1)
%             z_vals(i,j) = final_surface_roughness(end,node_counter);
%             node_counter = node_counter + 1;
%         end
%     end
    %=====================================================================
%     plotno = 2;
%     tp_true = 1;
%     plot_contour_tribo(plotno,X,Y,z_vals,num_levels,tp_true,tpos)
    %==========================================================================
    %==========================================================================
    i2 = find(i_tot_all, 1, 'last');
    i_tot = i_tot_all(1:i2);
    
%     i1 = find(t_all, 1, 'last');    
    t = t_all(1:i2);

    
    modifier_area = contact_area; %((surface_dimensions(1)*1.0e-7))*((surface_dimensions(2)*1.0e-7)); %cm2 %/nodes(1) /nodes(2)
    modifier_nodes = nodes(1)*nodes(2);
    modifier_A_to_fA = 1.0e15;

    % Create csv output files
    MatrixOutput1 = zeros(length(t), 2); 
    MatrixOutput1(:,1) = t;
    MatrixOutput1(:,2) = i_tot.*(modifier_area*modifier_A_to_fA);
%     
%     count = 1;
%     for i = 1:total_runs
%         for j = 1:number_of_time_steps
%             MatrixOutput1(j,count) = t(i,j);
% %             MatrixOutput1(j,count+1) = i_tot(i,j); %A/node * m^2
%             MatrixOutput1(j,count+1) = (i_tot(i,j)*modifier_area*modifier_A_to_fA); %fA  modifier_nodes* *modifier_area 
%         end        
%         count = count + 2;
%     end
       
    writematrix(MatrixOutput1,output_files,'Delimiter','comma') 
    %==========================================================================
    %=====================================================================
    % Plot formatting
    tick_label_size = 16;
    axis_label_size = 18;
    title_label_size = 20;
    plot_line_width = 3;
    axis_line_width = 2;
    marker_size = 3;
    font_weight = 'bold';
    color_range = {'b','b','b','b','b','b','b','b','b', ...
        'r','r','r','r','r','r','r','r','r', ...
        'g','g','g','g','g', ...
        'k','k','k','k','k', ...
        'c'};
    marker_vals_data = {'s','o','+','s','o','+','s','o','+','s','o','+', ...
        's','o','+','s','o','+','s','o','+'};
    % =====================================================================       
    %==========================================================================    
    data_time = data(:,1);
    data_current = data(:,2);
    
    if fig_num == 2
        t = t + 10.0;
    elseif fig_num == 3
        t = t + 8.5;
    end
    
    figure(fig_num)
    hold on
    axis square    

    mv = marker_vals_data{fig_num};
    color_vals = color_range{fig_num};        
    plot(t, i_tot.*(modifier_area*modifier_A_to_fA), strcat(color_vals,mv), 'MarkerSize',marker_size,'LineWidth',plot_line_width) %,mv '-', .*modifier_area .*modifier_nodes 
    
    color_vals2 = color_range(fig_num+9);
    plot(data_time, data_current, '-r') %, strcat('-',color_vals2), 'MarkerSize',marker_size,'LineWidth',plot_line_width-1
    
    xlabel('t (s)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    ylabel('i (fA)', 'FontSize', axis_label_size,'FontWeight',font_weight) 
    
    box on
    ax = gca;

    ax.FontSize = tick_label_size;
    ax.FontWeight = font_weight;
    ax.LineWidth = axis_line_width;
    ax.XMinorTick = 'on';    
    hold off    
    
    fig_name = strcat('Figure ',num2str(fig_num));
    saveas(fig_num,fig_name,'png');
    
end

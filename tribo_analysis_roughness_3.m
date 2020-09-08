function tribo_analysis_roughness_3(zsim)
    %==========================================================================
    % This function accetps a simulation structure, plots the initial and
    % final morphologies, if required, and calls the AFM_scan function with
    % the necessary parameters.
    %==========================================================================

    %==========================================================================
    % Simulation of tribocorrosion damage
    %========================================================================== 
    data = csvread(char(zsim.exp_filename));

    x_pos = 0.0:zsim.delta_xyz(1):zsim.oxide_xyz(1);
    y_pos = 0.0:zsim.delta_xyz(2):zsim.oxide_xyz(2); 
    
    %=====================================================================
    % Set-up initial surface morphology and plot it
    %=====================================================================    
    initial_surface_roughness = build_tribo_surface(zsim.morphology, zsim.N_asp, x_pos, y_pos, zsim.oxide_xyz(3), zsim.oxide_xyz(4));
    
    maz = max(initial_surface_roughness);
    miz = min(initial_surface_roughness);
    peak_height = max(maz);
    peak_valley = min(miz);
    delta_height = (peak_height - peak_valley)/15;
    
    % Set up plotting the initial surface morphology
    show_init_morph = 0;
    if show_init_morph == 1
        num_levels = peak_valley:delta_height:peak_height;    
        plotno = 1;
        tp_true = 0;
        tpos = 0.0;
        [X,Y] = meshgrid(y_pos, x_pos);
        plot_contour_tribo(plotno,X,Y,initial_surface_roughness,num_levels,tp_true,tpos)        
    end

    %=====================================================================
    
    %=====================================================================
    % Call AFM_scan_3 to model the current from the AFM tip interaction
    % with the oxide
    %===================================================================== 
    % Adjust dt to match the multiplier from multi_AFM_scan
    actual_dt = zsim.delta_xyz(1)/zsim.v_tip; 
    
    number_passes = 1;
    number_of_intervals = number_passes*round(zsim.scanning_time/actual_dt); %round(nodes(1) * nodes(2)); %(round (1.6*(nodes(1) * nodes(2))) + 10); %        
    relaxation_time_intervals = round(zsim.relaxation/actual_dt); 
    number_of_time_steps = number_of_intervals+relaxation_time_intervals+1;

    i0p2 = zsim.i0_pass/(zsim.comp_nodes(1) * zsim.comp_nodes(2));
    
    % Perform the AFM tip scan
    [tpos, final_surface_roughness, i_tot_all, t_all, contact_area] = ...
        AFM_scan_3(x_pos, y_pos, initial_surface_roughness, zsim.comp_nodes, ...
        number_of_intervals, relaxation_time_intervals, number_of_time_steps, actual_dt, ...
        'Diamond', 'Cr2O3', zsim.load, zsim.v_tip, zsim.eta, zsim.alpha, ...
        zsim.i0_field, zsim.Ef, zsim.i0_mono, i0p2);
    %=====================================================================
    
    %=====================================================================
    % Update the morphology of the asperities for plotting
    %===================================================================== 
    show_final_morph = 0;
    if show_final_morph == 1
        node_counter = 1;
        z_vals = zeros(size(initial_surface_roughness));
        for i = 1:(nodes(1)+1)
            for j = 1:(nodes(2)+1)
                z_vals(i,j) = final_surface_roughness(end,node_counter);
                node_counter = node_counter + 1;
            end
        end
        %=====================================================================
        plotno = 2;
        tp_true = 1;
        plot_contour_tribo(plotno,X,Y,z_vals,num_levels,tp_true,tpos)
        %=====================================================================        
    end
    %=====================================================================
    % Create CSV output files for use with other plotting software
    %=====================================================================
    % Remove any trailing zeros from the output arrays from AFM_scan_3
    i2 = find(i_tot_all, 1, 'last');
    i_tot = i_tot_all(1:i2); 
    t = t_all(1:i2);
    
    modifier_area = contact_area; %((surface_dimensions(1)*1.0e-7))*((surface_dimensions(2)*1.0e-7)); %cm2 %/nodes(1) /nodes(2)
%     modifier_nodes = zsim.comp_nodes(1)*zsim.comp_nodes(2);
    modifier_A_to_fA = 1.0e15;
    
    % Setup the output matrices
    MatrixOutput1 = zeros(length(t), 2); 
    MatrixOutput1(:,1) = t;
    MatrixOutput1(:,2) = i_tot.*(modifier_area*modifier_A_to_fA); 
    % Create the output file and write the output matrix to it
    writematrix(MatrixOutput1,char(zsim.output_filename),'Delimiter','comma') 
    %=====================================================================
    %=====================================================================
    % Plot the model current vs the measured current for comparison
    % purposes
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
    %=====================================================================     
    %=====================================================================    
    data_time = data(:,1);
    data_current = data(:,2);
    
    if zsim.figure_number == 2
        t = t + 10.0;
    elseif zsim.figure_number == 3
        t = t + 8.5;
    end
    
    figure(zsim.figure_number)
    hold on
    axis square    

    mv = marker_vals_data{zsim.figure_number};
    color_vals = color_range{zsim.figure_number};        
    plot(t, i_tot.*(modifier_area*modifier_A_to_fA), strcat(color_vals,mv), 'MarkerSize',marker_size,'LineWidth',plot_line_width) %,mv '-', .*modifier_area .*modifier_nodes 
    
    color_vals2 = color_range(zsim.figure_number+9);
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
    
    %=====================================================================
    % Save the figure as a PNG file
    %=====================================================================
    fig_name = strcat('Figure ',num2str(zsim.figure_number));
    saveas(zsim.figure_number,fig_name,'png');
    %=====================================================================
    
end

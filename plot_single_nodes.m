function plot_single_nodes
    clc;
    clear all;
    %=====================================================================
    % Plot formatting
    tick_label_size = 16;
    axis_label_size = 18;
    title_label_size = 20;
    plot_line_width = 3;
    axis_line_width = 2;
    marker_size = 7;
    font_weight = 'bold';
    color_range = {'b','r','g','k','c', ...
        'b','r','g','k','c', ...
        'b','r','g','k','c', ...
        'b','r','g','k','c', ...
        'b'};
    marker_vals_data = {'s','o','+','s','o','+','s','o','+','s','o','+', ...
        's','o','+','s','o','+','s','o','+'};
    %=====================================================================
    
    %=====================================================================
    % Names of the CSV output files. Will be stored in the same folder as
    % the execution code
    %=====================================================================
    directory = 'C:\Users\steve\OneDrive\Tribcorrosion\Matlab Scripts\Tribo-Model\SingleNode\';
    output_files = {'Model58GPa.csv', 'Model86GPa.csv', ...
        'Model103GPa.csv', 'Model1125GPa.csv', 'Model1259GPa.csv'};    
    %=====================================================================
    %=====================================================================
    nfiles = length(output_files);
    figure(30)
    hold on
    for i = 1:nfiles
        of = strcat(directory,output_files(i));
        mdl_res = csvread(char(of));
        time = mdl_res(:,1);
        I_vals = mdl_res(:,2);    
        color_vals = color_range{i};
        marker_vals = marker_vals_data(i);
        ptype = char(strcat('-',color_vals)); %,marker_vals
        plot(time,I_vals,ptype, 'MarkerSize',marker_size,'LineWidth',plot_line_width)
    end
    
    xlabel('Time (s)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    ylabel('I (fA)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    
    xlim([0.0 100.0])
    
    axis square    
    box on
    ax = gca;
    ax.FontSize = tick_label_size;
    ax.FontWeight = font_weight;
    ax.LineWidth = axis_line_width;
    ax.XMinorTick = 'on';   
    
    hold off
    
end
function plot_contour_tribo(plotno,X,Y,Z,num_levels,tp_true,tpos)
    %=====================================================================
    % Plot formatting
    tick_label_size = 16;
    axis_label_size = 18;
    title_label_size = 20;
    plot_line_width = 3;
    axis_line_width = 2;
    marker_size = 8;
    font_weight = 'bold';
    color_range = {'b','b','b','b','b','b','b','b','b', ...
        'r','r','r','r','r','r','r','r','r', ...
        'g','g','g','g','g', ...
        'k','k','k','k','k', ...
        'c'};
    marker_vals_data = {'s','o','+','s','o','+','s','o','+','s','o','+', ...
        's','o','+','s','o','+','s','o','+'};
    % =====================================================================     
    
    figure(plotno)
    hold on
    axis equal
    contour(X,Y,Z,num_levels,'LineWidth',plot_line_width)
    
    if tp_true == 1
        plot(tpos(:,1),tpos(:,2),'r+','MarkerSize',marker_size,'LineWidth',plot_line_width-1)
        plot(tpos(:,1),tpos(:,3),':r','LineWidth',plot_line_width-2)
        plot(tpos(:,1),tpos(:,4),':r','LineWidth',plot_line_width-2)        
    end
    
%     zlim([0 50])
    stopx = X(end,end);
    stopy = Y(end,end);
    xlim([0 stopx])
    ylim([0 stopy])
    
    xlabel('y (nm)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    ylabel('x (nm)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    zlabel('z (nm)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    
    box on
    ax = gca;

    ax.FontSize = tick_label_size;
    ax.FontWeight = font_weight;
    ax.LineWidth = axis_line_width;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on'; 
    
    c = colorbar;
    c.Label.String = 'Elevation (nm)';    
    
%     view(0,90)
    hold off   
  
end
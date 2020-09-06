function film_growth
    clc;
    clear all;
    
    Faraday_Constant = 96485; %coul/mol
    R = 8.314; %J/mol K
    T = 298; %K
    Eapp = [0.2, 0.4]; %VSCE
    Ecorr = -0.3; %VSCE
    eta = Eapp - Ecorr;
    
%     alpha = 0.25; %0.68; %0.25; %0.5;
%     i0_growth = 1.0e-16; %A
%     k_film = 1.0e10; %5.0e7; %um/A.s

    alpha = 0.1; %0.16; %0.25; %0.25; %0.68; %0.25; %0.5;
    i0_growth = 1.0e-9; %1.0e-2 * 1.0e-18; %1.0e-13; %1.0e-16; %A
    k_film = 5.0e7; %1.0e9; %1.0e10; %5.0e7; %um/A.s  
    
    z = 2;
    g_plus = (alpha * z * Faraday_Constant)/ ( R * T);
    
    E_film = Eapp(1)/(2.0e-3); %20; %100; %500; %250; %V/um    
    
    t = 0.0:0.01:1.0;

    num = i0_growth * exp(g_plus.*eta(1));
    denom = 1.0 + ((g_plus * E_film *k_film * num) .* t);
    
    damage_current = (num./denom);    
    
    %=====================================================================
    % Plot formatting
    tick_label_size = 16;
    axis_label_size = 18;
    title_label_size = 20;
    plot_line_width = 3;
    axis_line_width = 2;
    marker_size = 8;
    font_weight = 'bold';
    marker_vals_data = {'s','o','+'};
    % =====================================================================  
    
    figure(3)
    hold on
    plot(t, damage_current, '-r','LineWidth',plot_line_width)
    
    xlabel('t (s)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    ylabel('i_{damage} (A)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    
    ylim([1.0e-12 1.0e-4])
    
    axis square
    box on 
    
    ax = gca;

    ax.FontSize = tick_label_size;
    ax.FontWeight = font_weight;
    ax.LineWidth = axis_line_width;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.YScale = 'log';
    
    hold off
end
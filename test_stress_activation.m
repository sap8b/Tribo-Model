function test_stress_activation
    clc; 
    clear all;
    %=====================================================================
    % Experimental data files to be loaded - must be in same folder as the
    % execution code
    %=====================================================================
    exp_files = {'Baseline58GPa.csv', 'Baseline86GPa.csv', ...
        'Baseline103GPa.csv', 'Baseline1125GPa.csv', 'Baseline1295GPa.csv'};
    %=====================================================================
    % Names of the CSV output files. Will be stored in the same folder as
    % the execution code
    %=====================================================================
    output_files = {'Model58GPa.csv', 'Model86GPa.csv', ...
        'Model103GPa.csv', 'Model1125GPa.csv', 'Model1259.csv'};
    %=====================================================================
    nfiles = length(exp_files);
    I_avg = zeros(1,nfiles);
    I_avg_with_cutoff = zeros(1,nfiles);
    for i = 1:nfiles
        data = csvread(char(exp_files(i)));
        I_avg(1,i) = mean(data(:,2));
        mdl_res = csvread(char(output_files(i)));
        I_avg_with_cutoff(1,i) = mean(mdl_res(:,2));
    end
    %=====================================================================
    % Plot formatting
    tick_label_size = 16;
    axis_label_size = 18;
    title_label_size = 20;
    plot_line_width = 3;
    axis_line_width = 2;
    marker_size = 7;
    font_weight = 'bold';
    color_range = {'b','b','b','b','b','b','b','b','b', ...
        'r','r','r','r','r','r','r','r','r', ...
        'g','g','g','g','g', ...
        'k','k','k','k','k', ...
        'c'};
    marker_vals_data = {'s','o','+','s','o','+','s','o','+','s','o','+', ...
        's','o','+','s','o','+','s','o','+'};
    %===================================================================== 
    
    stress = [5.8,  8.6, 10.3, 11.25, 12.59].*1.0e9; %GPa %
    stress_no86 = [5.8,  10.3, 11.25, 12.59].*1.0e9; %GPa 
%     I_avg = [33.03812, 54.12765, 56.66675, 80.4581, 122.44677].*1.0e-15; %A Experimental Data
    I_avg = I_avg.*1.0e-15; %A Experimental Data
    
    I_avg_model_nocutoff = [32.4623,  64.77407, 73.97724, 94.48909].*1.0e-15; %A Tribo Model - no cutoff
    
%     I_avg_with_cutoff = [29.86840232,  37.21755627, 47.79358577, 60.87252707, 71.82683551].*1.0e-15; %A Tribo Model - t_cutoff = 32s
    I_avg_with_cutoff = I_avg_with_cutoff.*1.0e-15; %A Tribo Model - t_cutoff = 32s
    
    fn = @exptest2;
    x1 = 5.0:0.1:13.0;
    x = x1.*1.0e9;    
    
    % Fit the exp data
    b0 = [14.0e-15, 1.0e-10 -0.1];    
    sp = [stress(3), stress(4), stress(5)];
    ip = [I_avg(3),I_avg(4), I_avg(5)];    
    mdl_tbc_600_data = fitnlm(sp,ip,fn,b0);
    
    disp(mdl_tbc_600_data)
    y = feval(mdl_tbc_600_data,x);
    
    b1 = [13.0e-15 1.0e-10 -0.1];
    mdl_tbc_600_no_cutoff = fitnlm(stress_no86,I_avg_model_nocutoff,fn,b1);
    
    disp(mdl_tbc_600_no_cutoff)    
    y2 = feval(mdl_tbc_600_no_cutoff,x);
    
    b2 = [15.0e-15 1.0e-10 -0.1];
    sp = [ stress(2), stress(3), stress(4), stress(5)];
    ip = [I_avg_with_cutoff(2), I_avg_with_cutoff(3), I_avg_with_cutoff(4), I_avg_with_cutoff(5)];
    mdl_tbc_600_cutoff32 = fitnlm(sp,ip,fn,b2);
    
    disp(mdl_tbc_600_cutoff32)    
    y3 = feval(mdl_tbc_600_cutoff32,x);  
    
    kb = 1.38e-23; %
    T = 298;
    Na = 6.02e23;
    conv_m3_A3 = 1.0e-10*1.0e-10*1.0e-10;
    
    exp_fit_Ea_coeff = mdl_tbc_600_data.Coefficients.Estimate(3);
    exp_fit_Va_coeff = mdl_tbc_600_data.Coefficients.Estimate(2);
    
    mdl_nocutoff_fit_Ea_coeff = mdl_tbc_600_no_cutoff.Coefficients.Estimate(3);
    mdl_nocutoff_fit_Va_coeff = mdl_tbc_600_no_cutoff.Coefficients.Estimate(2);
    
    mdl_cutoff_fit_Ea_coeff = mdl_tbc_600_cutoff32.Coefficients.Estimate(3);
    mdl_cutoff_fit_Va_coeff = mdl_tbc_600_cutoff32.Coefficients.Estimate(2);
    
    exp_fits = [(-exp_fit_Ea_coeff*Na)/1000, exp_fit_Va_coeff/conv_m3_A3].*(kb*T); %kJ/mol  A3
    mdl_nocutoff_fits = [(-mdl_nocutoff_fit_Ea_coeff*Na)/1000, mdl_nocutoff_fit_Va_coeff/conv_m3_A3].*(kb*T); %kJ/mol  A3
    mdl_cutoff_fits = [(-mdl_cutoff_fit_Ea_coeff*Na)/1000, mdl_cutoff_fit_Va_coeff/conv_m3_A3].*(kb*T); %kJ/mol  A3
    
    i_exp_str = strcat('i_{fit, Exp}, E_{a} = ', num2str(exp_fits(1)), ' kJ/mol; V_{a} = ', num2str(exp_fits(2)), 'A^{3}');
    i_mdl1_str = strcat('i_{fit, Model 1}, E_{a} = ', num2str(mdl_nocutoff_fits(1)), ' kJ/mol; V_{a} = ', num2str(mdl_nocutoff_fits(2)), 'A^{3}');
    i_mdl2_str = strcat('i_{fit, Model 2}, E_{a} = ', num2str(mdl_cutoff_fits(1)), ' kJ/mol; V_{a} = ', num2str(mdl_cutoff_fits(2)), 'A^{3}');    
    
    figure(1)
    hold on
    plot(stress,I_avg,'k+', 'MarkerSize',marker_size+2,'LineWidth',plot_line_width)
    plot(stress_no86,I_avg_model_nocutoff,'r^', 'MarkerSize',marker_size,'LineWidth',plot_line_width)
    plot(stress,I_avg_with_cutoff,'bo', 'MarkerSize',marker_size,'LineWidth',plot_line_width)
    
    plot(x,y,'-k', 'MarkerSize',marker_size,'LineWidth',plot_line_width)
    plot(x,y2,'-r', 'MarkerSize',marker_size,'LineWidth',plot_line_width)
    plot(x,y3,'-b', 'MarkerSize',marker_size,'LineWidth',plot_line_width)
    
    axis square
    legend('i_{avg}, Measured','i_{avg}, Model 1 - no cutoff', 'i_{avg}, Model 2 - \tau_{cutoff} = 32 s', ...
        i_exp_str,i_mdl1_str,i_mdl2_str,'Location','northwest')
    
    xlabel('Stress (Pa)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    ylabel('I (A)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    
%     title('i_{0,passive}', 'FontSize', title_label_size,'FontWeight',font_weight)
    
    box on
    ax = gca;
    ax.FontSize = tick_label_size;
    ax.FontWeight = font_weight;
    ax.LineWidth = axis_line_width;
    ax.XMinorTick = 'on';  
    legend boxoff
    hold off    
    
    hold off
    
    
end

function y = exptest(b,x)
    kb = 1.38e-23;
    T = 298; %K
    den = kb * T;
    exp_term = b(3).*x;
    exp_term4 = b(2) + exp_term;
    exp_term2 = exp_term4./den;
    exp_term3 = exp(exp_term2); 
    y = b(1) * exp_term3;
end

function y = exptest2(b,x)
    b(1) = 13.2e-15; %A
%     b(3) = -0.14886; %-2.0196595; %-0.242;
    y = b(1).*exp((b(2).*x) + b(3));
end

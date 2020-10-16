function test_stress_activation
    clc; 
    clear all;
    
    dir = 'C:\Users\steve\OneDrive\Tribcorrosion\Matlab Scripts\Tribo-Model\'; %HoldOutputs\
    kb = 1.38e-23; %
    T = 298;
    Na = 6.02e23;
    conv_m3_A3 = 1.0e-10*1.0e-10*1.0e-10;
    v_act_test = [(0.277*conv_m3_A3),(0.277*conv_m3_A3),(0.261*conv_m3_A3)];
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
        'Model103GPa.csv', 'Model1125GPa.csv', 'Model1259GPa.csv'};
    %=====================================================================
    nfiles = length(exp_files);
    I_avg = zeros(1,nfiles);
    I_avg_with_cutoff = zeros(1,nfiles);
    for i = 1:nfiles
        data = csvread(strcat(dir,char(exp_files(i))));
        I_avg(1,i) = mean(data(:,2));
        
        mdl_res = csvread(strcat(dir,char(output_files(i))));
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
%     I_avg = I_avg.*1.0e-15; %A Experimental Data
    I_avg = [37.5, 89.75, 86.6, 104.57, 155.18].*1.0e-15; %A Experimental Data
    
    I_avg_model_nocutoff = [32.4623,  64.77407, 73.97724, 94.48909].*1.0e-15; %A Tribo Model - no cutoff
    
    I_avg_with_cutoff_top = [29.86840232,  100.2, 114.8, 124.5, 139.9].*1.0e-15; %A Tribo Model - t_cutoff = 32s
    I_avg_with_cutoff = I_avg_with_cutoff.*1.0e-15; %A Tribo Model - t_cutoff = 32s
%     I_avg_with_cutoff = [33.82, 100.71, 115.37, 125.19, 140.78].*1.0e-15; %A Tribo Model - t_cutoff = 48s
        
    fn = @exptest2;
    x1 = 5.0:0.1:13.0;
    x = x1.*1.0e9;    
    
    %=====================================================================
    % Fit the exp data
    %=====================================================================
%     b0 = [10.0e-15, v_act_test(1) 0.3];    
%     sp = [stress(1), stress(2), stress(3), stress(4)];
%     ip = [I_avg(1),I_avg(2), I_avg(3), I_avg(4)];    
%     mdl_tbc_600_data = fitnlm(sp,ip,fn,b0);
%     
%     disp(mdl_tbc_600_data)
%     y = feval(mdl_tbc_600_data,x);
    %=====================================================================
    %=====================================================================
    % Fit the cutoff top model
%     b1 = [1.0e-14 9.0e-31 ((0.01/Na)*1000)/(kb*T)]; %[14.485e-15 v_act_test(2) -0.1];
%     sp = [stress(2), stress(3), stress(4), stress(5)]; %
%     ip = [I_avg_with_cutoff_top(2), I_avg_with_cutoff_top(3), I_avg_with_cutoff_top(4), I_avg_with_cutoff_top(5)]; %    
%     mdl_tbc_600_cutoff_top = fitnlm(sp,ip,fn,b1);
%     
%     disp(mdl_tbc_600_cutoff_top)    
%     y2 = feval(mdl_tbc_600_cutoff_top,x);
    %=====================================================================
    %=====================================================================
    % Fit the cutoff model    
    %=====================================================================
%     b2 = [1.0e-15 1.4474e-30 0.00040393];
%     sp = [ stress(1), stress(2), stress(3), stress(4), stress(5)]; %
%     ip = [I_avg_with_cutoff(1), I_avg_with_cutoff(2), I_avg_with_cutoff(3), I_avg_with_cutoff(4), I_avg_with_cutoff(5)]; %
%     mdl_tbc_600_cutoff32 = fitnlm(sp,ip,fn,b2);
%     
%     disp(mdl_tbc_600_cutoff32)    
%     y3 = feval(mdl_tbc_600_cutoff32,x);  
    %=====================================================================
    %=====================================================================
    % Linear fit
%     p = polyfit(sp,ip,1);
%     y4 = polyval(p,x);
    %=====================================================================
        
%     exp_fit_Ea_coeff = mdl_tbc_600_data.Coefficients.Estimate(3);
%     exp_fit_Va_coeff = mdl_tbc_600_data.Coefficients.Estimate(2);
% %     
% %     mdl_nocutoff_fit_Ea_coeff = mdl_tbc_600_cutoff_top.Coefficients.Estimate(3);
% %     mdl_nocutoff_fit_Va_coeff = mdl_tbc_600_cutoff_top.Coefficients.Estimate(2);
%     
%     mdl_cutoff_fit_Ea_coeff = mdl_tbc_600_cutoff32.Coefficients.Estimate(3);
%     mdl_cutoff_fit_Va_coeff = mdl_tbc_600_cutoff32.Coefficients.Estimate(2);
%     
%     exp_fits = [(-exp_fit_Ea_coeff*Na)/1000, (exp_fit_Va_coeff/conv_m3_A3)/(kb*T)].*(kb*T); %kJ/mol  A3
% %     mdl_nocutoff_fits = [(-mdl_nocutoff_fit_Ea_coeff*Na)/1000, mdl_nocutoff_fit_Va_coeff/conv_m3_A3].*(kb*T); %kJ/mol  A3
%     v_act_fit = (mdl_cutoff_fit_Va_coeff/conv_m3_A3)/(kb*T);
%     mdl_cutoff_fits = [(-mdl_cutoff_fit_Ea_coeff*Na)/1000, v_act_fit].*(kb*T); %kJ/mol  A3
    
    i0_mdl = 1.0e-14;
    v_act_mdl = 9.0e-31;
    e_act_act = ((0.01/Na)*1000)/(kb*T);
    s_term = (x.*v_act_mdl)./(kb*T); 
    e_term = (-e_act_act + s_term);
    expval = exp(e_term);
    i_mdl_guess_red = i0_mdl.*expval;
    
    i0_mdl = 2.1434e-15; %1.0e-15;
%     v_act_mdl = 17.0e-31; %1.822e-30;
%     e_act_act = 1.0e-4; %0.00040393; %((0.01/Na)*1000)/(kb*T);
%     s_term = (x.*v_act_mdl)./(kb*T); 
%     e_term = (-e_act_act + s_term);
    e_term = 3.4310E-10.*x;
    expval = exp(e_term);
    i_mdl_guess_blue = i0_mdl.*expval;
    
    %=====================================================================
    % Setup the output matrices
    MatrixOutput1 = zeros(length(stress),4);
    MatrixOutput1(:,1) = stress;
    MatrixOutput1(:,2) = I_avg;
    MatrixOutput1(:,3) = stress;
    MatrixOutput1(:,4) = I_avg_with_cutoff;
    writematrix(MatrixOutput1,'BestFitExpMdl.xls','Sheet',1)
    
    MatrixOutput2 = zeros(length(x), 2); 
    MatrixOutput2(:,1) = x;
    MatrixOutput2(:,2) = i_mdl_guess_red; 
    % Create the output file and write the output matrix to it
    writematrix(MatrixOutput2,'BestFitExpMdl.xls','Sheet',2)    
    %=====================================================================
    
%     i_exp_str = strcat('i_{fit, Exp}, E_{a} = ', num2str(exp_fits(1)), ' kJ/mol; V_{a} = ', num2str(exp_fits(2)), 'A^{3}');
%     i_mdl1_str = strcat('i_{fit, Model 1}, E_{a} = ', num2str(mdl_nocutoff_fits(1)), ' kJ/mol; V_{a} = ', num2str(mdl_nocutoff_fits(2)), 'A^{3}');
%     i_mdl2_str = strcat('i_{fit, Model 2}, E_{a} = ', num2str(mdl_cutoff_fits(1)), ' kJ/mol; V_{a} = ', num2str(mdl_cutoff_fits(2)), 'A^{3}');    
    
    figure(20)
    hold on
    plot(stress,I_avg,'r+', 'MarkerSize',marker_size+2,'LineWidth',plot_line_width)
%     plot(stress_no86,I_avg_model_nocutoff,'r^', 'MarkerSize',marker_size,'LineWidth',plot_line_width)
    plot(stress,I_avg_with_cutoff,'bo', 'MarkerSize',marker_size,'LineWidth',plot_line_width)
%     plot(stress,I_avg_with_cutoff_top,'g^','MarkerSize',marker_size,'LineWidth',plot_line_width)
    
%     plot(x,y,':k', 'MarkerSize',marker_size,'LineWidth',plot_line_width)
%     plot(x,y2,'-g', 'MarkerSize',marker_size,'LineWidth',plot_line_width)
    
%     plot(x,y4,'-.b','MarkerSize',marker_size,'LineWidth',plot_line_width)
    
    plot(x,i_mdl_guess_red,'-r', 'MarkerSize',marker_size,'LineWidth',plot_line_width)
    plot(x,i_mdl_guess_blue,'-b', 'MarkerSize',marker_size,'LineWidth',plot_line_width)
    
    axis square
%     legend('i_{avg}, Measured','i_{avg}, Model 1 - no cutoff', 'i_{avg}, Model 2 - \tau_{cutoff} = 48 s', ...
%         i_exp_str,i_mdl1_str,i_mdl2_str,'Location','northwest')
    
%     legend('i_{avg}, Measured','i_{avg}, Model 2 - \tau_{cutoff} = 7500 s', ...
%         i_exp_str,i_mdl2_str, 'i_{linear}','i_{guess}','Location','northwest')
    
    xlabel('Stress (Pa)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    ylabel('I (A)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    
%     title('i_{0,passive}', 'FontSize', title_label_size,'FontWeight',font_weight)
    
    box on
    ax = gca;
    ax.FontSize = tick_label_size;
    ax.FontWeight = font_weight;
    ax.LineWidth = axis_line_width;
    ax.XMinorTick = 'on';  
%     legend boxoff
    
    hold off
    
    
end

function y = exptest(b,s)
    b(1) = 14.485e-15; %A
    b(3) = 0.1/(1.0e-10*1.0e-10*1.0e-10);
    kb = 1.38e-23;
    T = 298; %K
    den = kb * T;
    exp_term = b(3).*s;
    exp_term4 = b(2) + exp_term;
    exp_term2 = exp_term4./den;
    exp_term3 = exp(exp_term2); 
    y = b(1) * exp_term3;
end

function y = exptest2(b,x)
    i0 = b(1); %14.485e-15; %A
    
    kb = 1.38e-23; %
    T = 298;
    v_act = b(2); %(0.1*conv_m3_A3);
    stress_effect =  (v_act.* x)./(kb*T);

    E_act = b(3);

    y = i0.*exp(stress_effect + E_act);
end



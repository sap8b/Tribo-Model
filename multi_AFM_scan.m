function multi_AFM_scan
    %=====================================================================
    % The purpose of this model is to simulate the corrosion current
    % measured from the tribological contact between a non-conductive AFM
    % tip and a polarized metal surface immmersed in a corrosive
    % electrolyte.
    % 
    % This function serves as the top-level function that defines the
    % experimental file names to load the experimental data for comparison,
    % the simulation output file names for saving simulation results, the
    % model parameters as arrays with individual entries corresponding to
    % the relevant simulation conditions, and the principal execution loop
    % that calls the tribo_analysis_roughness with the parameters needed
    %
    % Revision Record
    %
    % Version       Date                Programmer
    % -------       ----                ----------
    % 1.0           17 August 2020      S.A. Policastro
    % 2.0           06 September 2020   S.A. Policastro
    %    
    %=====================================================================
    clc;
    clear all;

    tic
    
    %=====================================================================
    % Order of scans: 5.8, 8.6, 10.3, 11.25, 12.95 GPa
    %=====================================================================
    total_exp_time = [86.0, 160.0, 160.0, 170.0, 185.0]; %s
    time_scanning =  [61.0, 110.0, 134.0, 142.0, 157.0]; %s
    relaxation_time = total_exp_time - time_scanning;
    %=====================================================================
    %=====================================================================
    % Experimental data files
    %=====================================================================
    exp_files = {'Baseline58GPa.csv', 'Baseline86GPa.csv', ...
        'Baseline103GPa.csv', 'Baseline1125GPa.csv', 'Baseline1295GPa.csv'};
    output_files = {'Model58GPa.csv', 'Model86GPa.csv', ...
        'Model103GPa.csv', 'Model1125.csv', 'Model1295.csv'};
    %=====================================================================
    %=====================================================================
    % Initialize variables
    %=====================================================================
    % Define probe surface
    % [x y zmin zmax]
    surface_dimensions = [20000, 15000, 2.6, 2.4]; % nm
    v_tip = 4000.0; % nm/s
    dt = 0.01; %s
    base_dx = v_tip * dt; %nm
    dx = 4*base_dx;
    % Create computational cell
    % [NX NY, NZ]
%     nodes = [80, 60, 10];
    nodes = [round(surface_dimensions(1)/dx), round(surface_dimensions(2)/dx), 10];
    % Create position deltas
    % [dx dy dz]
    delta_pos = [(surface_dimensions(1)+1)/nodes(1), (surface_dimensions(2)+1)/nodes(2), (surface_dimensions(3) - surface_dimensions(4))/nodes(3)]; %nm  
    
    node_area_nm2 = delta_pos(1) * delta_pos(2); %nm2
    node_area_cm2 = node_area_nm2 * 1.0e-7 * 1.0e-7;
    %=====================================================================
    % Set up the tip-to-oxide interaction parameters
    %=====================================================================
    
    L_base = [0.009, 0.03, 0.0495, 0.064, 0.09]; %
    L_applied = L_base.*1.0e-5; % N
    %=====================================================================
    %=====================================================================
    % Setup the electrochemical system parameters
    %=====================================================================    
    Eapp = 0.6; % 0.54; %[0.0, 0.2, 0.4, 0.54]; %VSCE
    Ecorr = -0.3; %VSCE
    eta = Eapp - Ecorr; %V
    %=====================================================================    
    %=====================================================================
    % Morphologies
    %=====================================================================
    surface_morphology = 1; %2; %1; %2; %2; %3; %4; 2; %3; %
    %=====================================================================
    %=====================================================================
    % Setup other parameters
%     alpha0 = [0.32, 0.385, 0.40, 0.42, 0.421];
%     alpha0 = [0.32, 0.40, 0.49, 0.49, 0.51];
    alpha0 = [0.33, 0.37, 0.37, 0.33, 0.33];

    i0_monolayer = [0.3, 1.0, 1.0, 0.3, 0.3] .*1.0e-3; %A/cm2
%     i0_growth_base = [12.0, 2.0, 0.07, 0.06, 0.03] .* 1.0e-10; %A/cm2
%     i0_growth_base = [14.0, 3.0, 0.4, 14.0, 14.0] .* 1.0e-10; %A/cm2
    i0_growth_base = [5.483, 1.176, 0.6, 14.0, 14.0] .* node_area_cm2; %A/cm2
    
    i0_pass = [1.0, 1.0, 1.0, 1.0, 1.0].*1.0e-2;%A/cm2
%     E_film = [0.38, 0.6, 0.76, 0.80, 0.90];
    E_film = [0.8, 0.8, 0.8, 0.8, 0.8]; % V/nm
    %=====================================================================    
    number_of_scans = length(total_exp_time);
    
    for idx = 3:3 %number_of_scans
        tribo_analysis_roughness_3(total_exp_time(idx), time_scanning(idx), ...
            relaxation_time(idx), exp_files{idx}, output_files{idx}, ...
            L_applied(idx), eta, surface_morphology, alpha0(idx), ...
            i0_growth_base(idx), E_film(idx),i0_monolayer(idx), i0_pass(idx), ...
            surface_dimensions, delta_pos, nodes, v_tip, idx);
    end
    
    toc
end
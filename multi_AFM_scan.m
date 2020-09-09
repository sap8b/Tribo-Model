    %%
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
    % that calls the tribo_analysis_roughness_3 with the parameters needed.
    % 
    % The overall governing equations for the measured current are:
    %
    % i_{monolayer} - the current arising from the formation of the initial 
    % monolayer of the oxide following the AFM tip passage
    %
    % i_{growth} - the current arising from the high field film growth 
    % following the initial monolayer formation
    %
    % i_{passive} - the long-term background passive current once the oxide 
    % film has stabilized 
    %%
    % $i_{tot} = i_{monolayer} + i_{growth} + i_{passive}$
    %
    % $i_{monolayer} = i_{0,m}e^{-Kt}$
    %
    % $i_{growth} = i_{0,g}e^{\eta g^{+}}(1 +
    % k_{film} E_{film} g^{+} i_{0}e^{\eta g^{+}} \Delta t)^{-1}$
    %
    % $i_{passive} = i_{0,p}$   
    %%
    %=====================================================================
function multi_AFM_scan
    %=====================================================================
    % Revision Record
    %
    % Version       Date                Programmer
    % -------       ----                ----------
    % 1.0           17 August 2020      S.A. Policastro
    % 2.0           06 September 2020   S.A. Policastro
    % 3.0           08 September 2020   S.A. Policastro   
    %=====================================================================

    clc;
    clear all;
    
    % Start check for wall-clock time
    tic
    
    %=====================================================================
    % Order of scans: $\sigma_{max}$ = 5.8, 8.6, 10.3, 11.25, 12.95 GPa    
    %=====================================================================
    total_exp_time = [86.0, 160.0, 160.0, 170.0, 185.0]; %s
    time_scanning =  [61.0, 110.0, 132.0, 142.0, 157.0]; %s
    relaxation_time = total_exp_time - time_scanning;
    %=====================================================================
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
        'Model103GPa.csv', 'Model1125.csv', 'Model1295.csv'};
    %=====================================================================
    %=====================================================================
    % Initialize variables
    %=====================================================================
    % Define dimensions of the surface
    % [x y zmin zmax]
    surface_dimensions = [20000, 15000, 2.6, 2.4]; % nm
    % Define the AFM tip scanning velocity    
    v_tip = 4000.0; % nm/s
    % Define the time-step
    dt = 0.01; %s
    % Define the basic unit of length
    base_dx = v_tip * dt; %nm
    % Define a length multiplier to reduce the size of the surface scanning
    % loops and speed up code execution
    length_multiplier = 4;
    dx = length_multiplier*base_dx;
    % Create computational cell
    % [NX NY, NZ]
    nodes = [round(surface_dimensions(1)/dx), round(surface_dimensions(2)/dx), 10];
    % Create position deltas
    % [dx dy dz]
    delta_pos = [(surface_dimensions(1)+1)/nodes(1), (surface_dimensions(2)+1)/nodes(2), (surface_dimensions(3) - surface_dimensions(4))/nodes(3)]; %nm  
    % Calculate node areas
    node_area_nm2 = delta_pos(1) * delta_pos(2); %nm2
    node_area_cm2 = node_area_nm2 * 1.0e-7 * 1.0e-7; % cm2
    %=====================================================================
    % Define the AFM tip loads applied to the surface
    %=====================================================================    
    L_base = [0.09, 0.3, 0.495, 0.64, 0.9]; %uN
    L_applied = L_base.*1.0e-6; % N
    %=====================================================================
    %=====================================================================
    % Setup the electrochemical system parameters
    %=====================================================================    
    Eapp = [0.2, 0.6]; %VSCE
    Ecorr = -0.3; %VSCE
    % $\eta$ - represents the overpotential between the applied voltage and
    % the corrosion potential of the substrate
    eta = Eapp - Ecorr; %V
    %=====================================================================    
    %=====================================================================
    % Morphologies
    % 1 - Smooth surface with oxide thickness = 2.5 nm
    % 2 - Surface with "N_asperities" node-sized (point) asperities
    % 3 - Surface with 1 paraboloid asperity
    % 4 - Surface with "N_asperities" paraboloid asperities 
    %=====================================================================
    surface_morphology = 1; %2; %1; %2; %2; %3; %4; 2; %3; %
    N_asperities = 0;
    %=====================================================================
    %=====================================================================
    % Setup high-field film growth parameters for the different
    % experimetnal runs
    %=====================================================================   
    % $\alpsha$0 - represents the generalized transfer coefficient for the
    % rate-determining step of the oxide formation reaction
    %=====================================================================
%     alpha0 = [0.32, 0.385, 0.40, 0.42, 0.421];
%     alpha0 = [0.32, 0.40, 0.49, 0.49, 0.51];
%     alpha0 = [0.32, 0.37, 0.38, 0.392755, 0.405165];
    alpha0 = mdl_a_plus(L_base); 
    %=====================================================================
    % $i_{0,monolayer}$ - represents the current due to the formation of 
    % initial monolayer of oxide
    %=====================================================================
    i0_monolayer = [1.0, 1.0, 1.0, 1.0, 1.0] .*1.0e-3; %A/cm2
%     i0_growth_base = [12.0, 2.0, 0.07, 0.06, 0.03] .* 1.0e-10; %A/cm2
%     i0_growth_base = [14.0, 3.0, 0.4, 14.0, 14.0] .* 1.0e-10; %A/cm2
    %=====================================================================
    % $i_{growth}$ - represents the current from the high-field film growth
    %=====================================================================
%     i0_growth_base = [10, 1.28, 0.5, 0.3, 0.2] .* node_area_cm2; %A/cm2
    i0_growth_base = mdl_ig(L_base).* node_area_cm2; %A/cm2
    %=====================================================================
    % $i_{passive} - represents the current from the passivated surface at
    % long times
    %=====================================================================    
%     i0_pass = [1.0, 0.1, 0.01, 0.001976, 0.000103].*1.0e-2;%A/cm2
    i0_pass = mdl_ip(L_base); %A/cm2
    %=====================================================================
    % $E_{film}$ - represents the electric field that is set up across the
    % oxide layer between the metal surface and the electrolyte due to
    % sepration of positive metal ions and electrons
    %=====================================================================
%     E_film = [0.38, 0.6, 0.76, 0.80, 0.90];
%     E_film = [0.5, 0.8, 0.9, 0.93, 0.935]; % V/nm
    E_film = mdl_ef(L_base); % V/nm
    %=====================================================================
    %=====================================================================
    % Create a simulation parameter structure to assist with passing the
    % simulation parameters to the execution code
    %=====================================================================
    number_of_scans = length(total_exp_time);
    for i = 1:number_of_scans
        sim(i).total_time = total_exp_time(i);
        sim(i).scanning_time = time_scanning(i);
        sim(i).relaxation = relaxation_time(i);
        sim(i).exp_filename = exp_files(i);
        sim(i).output_filename = output_files(i);
        sim(i).load = L_applied(i);
        sim(i).eta = eta(2);
        sim(i).morphology = surface_morphology;
        sim(i).alpha = alpha0(i);
        sim(i).i0_field = i0_growth_base(i);
        sim(i).i0_mono = i0_monolayer(i);
        sim(i).i0_pass = i0_pass(i);
        sim(i).Ef = E_film(i);
        sim(i).oxide_xyz = surface_dimensions;
        sim(i).delta_xyz = delta_pos;
        sim(i).comp_nodes = nodes;
        sim(i).v_tip = v_tip;
        sim(i).N_asp = N_asperities;
        sim(i).figure_number = i;
    end
    
    
    for idx_iteration = 1:number_of_scans
        tribo_analysis_roughness_3(sim(idx_iteration));
    end
    %=====================================================================
    % Finish check for wall-clock time
    toc
end
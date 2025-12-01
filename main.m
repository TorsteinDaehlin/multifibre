%% Predictive Simulations of Human Gait

% This script starts the predictive simulation of human movement. The
% required inputs are necessary to start the simulations. Optional inputs,
% if left empty, will be taken from getDefaultSettings.m.

clear
close all
clc
% path to the repository folder
[pathRepo,~,~] = fileparts(mfilename('fullpath'));
% path to the folder that contains the repository folder
[pathRepoFolder,~,~] = fileparts(pathRepo);

%% Initialize S
addpath(fullfile(pathRepo,'DefaultSettings'))

[S] = initializeSettings('DHondt_et_al_2024_3seg');

%% Settings

% name of the subject
S.subject.name = 'DHondt_et_al_2024_3seg';

% path to folder where you want to store the results of the OCP
S.misc.save_folder  = fullfile(pathRepoFolder,'PredSimResults',[S.subject.name '_multifibre'],'tact_sensitivity', 'Walk');
% S.misc.save_folder  = fullfile(pathRepoFolder,'PredSimResults',[S.subject.name '_multifibre'],'fibre_shift','Run1_0_ms_new');
% S.misc.save_folder  = fullfile(pathRepoFolder,'PredSimResults',[S.subject.name], 'speeds');


% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
% S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
% S.solver.IG_selection = fullfile(S.misc.save_folder, 'DHondt_et_al_2024_3seg_v1.mot');

% S.solver.IG_selection = fullfile(pathRepoFolder,'PredSimResults',[S.subject.name '_multifibre'],'fibre_shift', ...
%     'Run_3_5_ms','QR_guess','DHondt_et_al_2024_3seg_56p_type_I.mot');
% S.solver.IG_selection_gaitCyclePercent = 200;
S.solver.IG_selection = 'quasi-random';
S.solver.IG_selection_gaitCyclePercent = 100;

% Set options for multi motor unit (MMU) muscle model
S.multifibre.use_multifibre_muscles = true;
S.multifibre.NFibre = 2;
S.multifibre.vMmax_range = [5 10]; % [6 10]; % Range of max contraction velocities as multiple of optimal fibre lengths
S.multifibre.tact_range = [0.025 0.045]; %[0.011 0.022]; % Range of activation time constants 
S.multifibre.beta = 0.6; % deactivation time constants are given by tact * (1 / beta).
S.multifibre.smeta = linspace(1.5, 2.5, S.multifibre.NFibre); % most efficient to least

% Set options for parameter shifts
fibre_type_shift = 1.0;
S.param_shift.slow_to_fast = fibre_type_shift;
% S.param_shift.fast_to_slow = fibre_type_shift;

% Set cost functional weights
S.weights.a = 1000; % Reduced to half of the original cost, as this weight is multiplied by a sum over twice as many activations.

% Set number of threads
S.solver.N_threads = 6;

S.misc.gaitmotion_type = 'HalfGaitCycle';

% Visualize bounds
S.misc.visualize_bounds = false;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

S.misc.forward_velocity = 1.33;

% S.metabolicE.model = 'MinettiAlexander';

% S.bounds.t_final.lower = 0.01;
S.bounds.default_coordinate_bounds = 'Running_Coordinate_Bounds.csv';
S.misc.task = 'walking';

% Run simulations as batch jobs, such that multiple simulations can run at
% the same time.
S.solver.run_as_batch_job = true;

% Set bounds on stride frequency (from Fukuchi et al. 2017 4.5 m/s [mean +/- SD])
% S.bounds.StrideFreq.lower = 1.306;
% S.bounds.StrideFreq.upper = 1.752;

%% Run predictive simulations
if S.solver.run_as_batch_job
    % fibre_type_shift = fliplr(linspace(0.1, 1, 5))';
    % fibre_type_shift = [0.775; 0.5500; 0.3250; 0.1000];
    % fibre_type_shift = 1;
    % n_meshes = [40 50 60 75 100 125];
    % speeds = [0.8 1.33 2.25 3.0 4.5]';
    tacts = [0.005 0.015; 0.015 0.025; 0.025 0.035; 0.025 0.045; 0.035 0.055; 0.045 0.065];
    % vmax = [0.5 0.75 1 1.25 1.5]';
    % vmax = [5 6 7 8 9]';
    % tacts = [0.05 0.015 0.025 0.035 0.045]';
    % tol = [2 3 4 5 6]';
    % w_Edot = [1 100 500 1000 2000]';
    for i = 1:size(tacts, 1)

        % S.misc.forward_velocity = speeds(i);
        % S.solver.N_meshes = n_meshes(i);
        % S.param_shift.slow_to_fast = fibre_type_shift(i);
        % S.param_shift.fast_to_slow = fibre_type_shift(i);
        S.multifibre.tact_range = tacts(i, :);
        % S.multifibre.vMmax_range = [5 10] .* vmax(i);
        % S.misc.custom_vMmax = vmax(i);
        % S.misc.custom_tact = tacts(i);
        % S.solver.tol_ipopt = tol(i);
        % S.weights.E = w_Edot(i);

        [savename] = runPredSim(S, osim_path);
        
    end
else
    [savename] = runPredSim(S, osim_path);
end

%% Plot results
% see .\PlotFigures\run_this_file_to_plot_figures.m for more

if ~S.solver.run_as_batch_job

    % set path to reference result
    % result_paths{1} = fullfile(pathRepo,'Tests','ReferenceResults',...
    %     'Falisse_et_al_2022','Falisse_et_al_2022_paper.mat');
    
    % set path to saved result
    result_paths{1} = fullfile(S.misc.save_folder,[savename '.mat']);
    
    % Cell array with legend name for each result
    legend_names = {'two_fibre'};
    
    % add path to subfolder with plotting functions
    addpath(fullfile(pathRepo,'PlotFigures'))
    
    figure_settings(1).name = 'Qs';
    figure_settings(1).dofs = {'all_coords'};
    figure_settings(1).variables = {'Qs'};
    figure_settings(1).savepath = [];
    figure_settings(1).filetype = {};

    % call plotting function
    plot_figures(result_paths, legend_names, figure_settings);

end

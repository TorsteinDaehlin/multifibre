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

[S] = initializeSettings('gait1018');

%% Settings

% name of the subject
S.subject.name = 'gait1018';

% path to folder where you want to store the results of the OCP
S.misc.save_folder  = fullfile(pathRepoFolder,'PredSimResults',[S.subject.name '_test']);

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.solver.IG_selection_gaitCyclePercent = 100;
% S.solver.IG_selection = 'quasi-random';

% Set options for multi motor unit (MMU) muscle model
S.multifibre.use_multifibre_muscles = true;
S.multifibre.NFibre = 2;
S.multifibre.vMmax_range = [6 10]; % Range of max contraction velocities as multiple of optimal fibre lengths
S.multifibre.tact_range = [0.01 0.02]; % Range of activation time constants 
S.multifibre.beta = 0.6; % deactivation time constants are given by tact * (1 / beta).
S.multifibre.smeta = linspace(1.5, 2.5, S.multifibre.NFibre); % most efficient to least

% Set cost functional weights
% S.weights.a = 1000; % Reduced to half of the original cost, as this weight is multiplied by a sum over twice as many activations.

% Set number of threads
S.solver.N_threads = 6;

S.misc.gaitmotion_type = 'HalfGaitCycle';

% Visualize bounds
% S.misc.visualize_bounds = true;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% Run simulations as batch jobs, such that multiple simulations can run at
% the same time.
S.solver.run_as_batch_job = false;

% S.misc.forward_velocity = 2.53;

% S.metabolicE.model = 'MinettiAlexander';

%% Run predictive simulations
if S.solver.run_as_batch_job
    % n_meshes = [40 50 60 75 100 125];
    % speeds = [0.83 1.33 1.83 2.33 2.58];
    for i = 1:length(speeds)

        S.misc.forward_velocity = speeds(i);
        % S.solver.N_meshes = n_meshes(i);

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
    legend_names = {'single fibre'};
    
    % add path to subfolder with plotting functions
    addpath(fullfile(pathRepo,'PlotFigures'))
    
    figure_settings(1).name = 'muscles_r';
    figure_settings(1).dofs = {'muscles_r'};
    figure_settings(1).variables = {'a'};
    figure_settings(1).savepath = [];
    figure_settings(1).filetype = {};

    % call plotting function
    plot_figures(result_paths, legend_names, figure_settings);

end

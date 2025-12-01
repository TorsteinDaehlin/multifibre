%% Predictive Simulations of Human Gait - Effect of triceps surae pennation

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
S.misc.save_folder  = fullfile(pathRepoFolder,'PredSimResults',['DHondt_et_al_2024_3seg_pennation'], 'scaled_lTs', 'gait_4_50', 'QR_guess');
S.misc.forward_velocity = 4.50;

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
% S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
% S.solver.IG_selection = fullfile(pathRepoFolder,'PredSimResults',['DHondt_et_al_2024_3seg_pennation'], 'gait_1_00', 'QR_guess','DHondt_et_al_2024_3seg_100_pennation.mot');
% S.solver.IG_selection_gaitCyclePercent = 200;
S.solver.IG_selection = 'quasi-random';
S.solver.IG_selection_gaitCyclePercent = 100;

% Set options for multi motor unit (MMU) muscle model
S.multifibre.use_multifibre_muscles = false;

% Set number of threads
S.solver.N_threads = 6;

S.misc.gaitmotion_type = 'HalfGaitCycle';

% Visualize bounds
S.misc.visualize_bounds = false;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% S.bounds.t_final.lower = 0.01;
S.bounds.default_coordinate_bounds = 'Running_Coordinate_Bounds.csv';
S.misc.task = 'walking';

% Run simulations as batch jobs, such that multiple simulations can run at
% the same time.
S.solver.run_as_batch_job = true;

% Load parameters
load(fullfile('C:\Users\nml-p\Documents\Projects\PredSimResults\DHondt_et_al_2024_3seg_pennation', ...
    'pennation_scale_factors.mat'), 'scales');

%% Run predictive simulations
if S.solver.run_as_batch_job
    for i = 1:length(scales.alphao)
        scale_C = {};
        for j = 1:length(scales.muscle_names)
            scale_C = [scale_C {scales.muscle_names(j), 'alphao', scales.alphao(i), ...
                scales.muscle_names(j), 'lTs', scales.lTs(j, i)}];
        end
        S.subject.scale_MT_params =  scale_C;
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

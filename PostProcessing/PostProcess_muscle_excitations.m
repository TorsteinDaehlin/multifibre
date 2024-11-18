function [R] = PostProcess_muscle_excitations(S, model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcess_muscle_excitations
%   This function computes muscle excitations from time derivative of muscle 
%   activations using Raasch's model.
%   More details in De Groote et al. (2009): DOI: 10.1080/10255840902788587
%
% INPUT:
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
%   - f_casadi -
%   * Struct containing all casadi functions.
%
%   - R -
%   * struct with simulation results
%
% OUTPUT:
%   - R -
%   * struct with simulation results
% 
% Original author: Lars D'Hondt
% Original date: 13/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


% Activation time constant
tact = reshape(ones(S.multifibre.NFibre,model_info.muscle_info.NMuscle).*model_info.muscle_info.tact', 1, S.multifibre.NFibre*model_info.muscle_info.NMuscle); 

% Deactivation time constant
tdeact = reshape(ones(S.multifibre.NFibre,model_info.muscle_info.NMuscle).*model_info.muscle_info.tdeact', 1, S.multifibre.NFibre*model_info.muscle_info.NMuscle); 

% Compute excitations
R.muscles.e = computeExcitationRaasch(R.muscles.a,R.muscles.da,tdeact,tact);


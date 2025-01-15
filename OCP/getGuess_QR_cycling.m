function guess = getGuess_QR_cycling(S,model_info,scaling,d)
% --------------------------------------------------------------------------
% getGuess_QR_cycling
%   This script provides an inital guess for the design variables for a cycling simulation.
%   The guess is quasi-random (QR). We set constant values to the muscle
%   variables and most joint variables. We only ensure
%   that the cycling speed is non=zero. We use a pre-defined final
%   time that is function of the imposed cycling frquency.
%   
% INPUT:
%   - S -
%   * setting structure S
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
%   - scaling -
%   * scale factors for all optimisation variables
% 
%   - d -
%   * degree of the interpolating polynomial of the collocation scheme
%
% OUTPUT:
%   - guess -
%   * initial guess values for all optimisation variables
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

N = S.solver.N_meshes; % number of mesh intervals
nq = model_info.ExtFunIO.jointi.nq;
NMuscle = model_info.muscle_info.NMuscle;
NFibre = S.multifibre.NFibre;
coordi = model_info.ExtFunIO.coordi;

%% Final time
% The final time is function of the imposed cycling frequency (given in
% rounds per minute). As we want the final time to be in seconds, we divide
% 60 s by the cycling frequency.
if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    guess.tf = 0.5 * (60 / S.cycling.rpm);
else
    guess.tf = 60 / S.cycling.rpm;
end

%% Qs
% The model is moving forward but with a standing position (Qs=0)
guess.Qs = zeros(N,nq.all);

%% Qdots
guess.Qdots = zeros(N,nq.all);

%% Qdotdots
guess.Qdotdots = zeros(N,nq.all);

%% Muscle variables
guess.a = 0.1*ones(N,NMuscle*NFibre);
guess.vA = 0.01*ones(N,NMuscle*NFibre);
guess.FTtilde = 0.1*ones(N,NMuscle);
guess.dFTtilde = 0.01*ones(N,NMuscle);
guess.SynH = 0.1*ones(N,NMuscle);
guess.SynW = 0.2;

%% Torque actuator activations
guess.a_a = 0.1*ones(N,nq.torqAct);
guess.e_a = 0.1*ones(N,nq.torqAct);

%% Pedal forces
n_FPedal = numel(model_info.ExtFunIO.input.Forces.pedal_force_r) + ...
    numel(model_info.ExtFunIO.input.Forces.pedal_force_l);
guess.FPedal = zeros(N, n_FPedal);

%% Add last mesh point to state variables
if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    % Lower limbs and trunk
    % Qs and Qdots are inverted after a half gait cycle BUT 
    % Pelvis: pelvis tilt and pelvis ty should be equal, pelvis
    % list, rot, tz should be opposite and pelvis tx should be equal
    % plus dist traveled.
    % Trunk: lumbar ext should be equal, lumbar bend and lumbar rot
    % should be of opposite.     
    % For "symmetric" joints, we invert right and left
    inv_X_Qs = zeros(1,nq.all);
    inv_X_Qs(1,model_info.ExtFunIO.symQs.QsInvA) = guess.Qs(1,model_info.ExtFunIO.symQs.QsInvB);
    inv_X_Qdots = zeros(1,nq.all);
    inv_X_Qdots(1,model_info.ExtFunIO.symQs.QdotsInvA) = guess.Qs(1,model_info.ExtFunIO.symQs.QdotsInvB);
    % For other joints, we take the opposite right and left
    inv_X_Qs(model_info.ExtFunIO.symQs.QsOpp) = ...
        -guess.Qs(1,model_info.ExtFunIO.symQs.QsOpp);           
    inv_X_Qdots(model_info.ExtFunIO.symQs.QsOpp) = ...
        -guess.Qdots(1,model_info.ExtFunIO.symQs.QsOpp);           
    guess.Qs = [guess.Qs; inv_X_Qs];
    guess.Qdots = [guess.Qdots; inv_X_Qdots];
    guess.a = [guess.a; guess.a(1,model_info.ExtFunIO.symQs.MActInvB)];
    guess.FTtilde = [guess.FTtilde; guess.FTtilde(1,model_info.ExtFunIO.symQs.MusInvB)];
    guess.a_a = [guess.a_a; guess.a_a(1,:)];
    guess.FPedal = [guess.FPedal; guess.FPedal(1,:)];
else
    guess.Qs = [guess.Qs; guess.Qs(1,:)];
    guess.Qdots = [guess.Qdots; guess.Qdots(1,:)];
    guess.a = [guess.a; guess.a(1,:)];
    guess.FTtilde = [guess.FTtilde; guess.FTtilde(1,:)];
    guess.a_a = [guess.a_a; guess.a_a(1,:)];
    guess.FPedal = [guess.FPedal; guess.FPedal(1,:)];
end


%% Scaling
guess.Qs = guess.Qs./repmat(scaling.Qs,N+1,1);
guess.Qdots = guess.Qdots./repmat(scaling.Qdots,N+1,1);
guess.Qdotdots  = guess.Qdotdots./repmat(scaling.Qdotdots,N,1);
guess.a         = (guess.a)./repmat(scaling.a,N+1,size(guess.a,2));
guess.FTtilde   = (guess.FTtilde)./repmat(scaling.FTtilde,N+1,1);
guess.vA        = (guess.vA)./repmat(scaling.vA,N,size(guess.vA,2));
guess.dFTtilde  = (guess.dFTtilde)./repmat(scaling.dFTtilde,N,size(guess.dFTtilde,2));
guess.FPedal    = (guess.FPedal)./repmat(scaling.FPedal, N+1, size(guess.FPedal, 2));
% guess.a_mtp_col = zeros(d*N,nq.mtp);
% guess.a_lumbar_col = zeros(d*N,nq.torso);

%% Collocation points
guess.a_col = zeros(d*N,NMuscle*NFibre);
guess.FTtilde_col = zeros(d*N,NMuscle);
guess.Qs_col = zeros(d*N,nq.all);
guess.Qdots_col = zeros(d*N,nq.all);
guess.a_a_col = zeros(d*N,nq.torqAct);
guess.dFTtilde_col = zeros(d*N,NMuscle);
guess.Qdotdots_col = zeros(d*N,nq.all);
guess.FPedal_col = zeros(d*N,n_FPedal);
for k=1:N
    guess.a_col((k-1)*d+1:k*d,:) = repmat(guess.a(k,:),d,1); 
    guess.FTtilde_col((k-1)*d+1:k*d,:) = repmat(guess.FTtilde(k,:),d,1);
    guess.Qs_col((k-1)*d+1:k*d,:) = repmat(guess.Qs(k,:),d,1);
    guess.Qdots_col((k-1)*d+1:k*d,:) = repmat(guess.Qdots(k,:),d,1);
    guess.a_a_col((k-1)*d+1:k*d,:) = repmat(guess.a_a(k,:),d,1);
    guess.dFTtilde_col((k-1)*d+1:k*d,:) = repmat(guess.dFTtilde(k,:),d,1);
    guess.Qdotdots_col((k-1)*d+1:k*d,:) = repmat(guess.Qdotdots(k,:),d,1);
    guess.FPedal_col((k-1)*d+1:k*d,:) = repmat(guess.FPedal(k,:),d,1);
end

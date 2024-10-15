function test_multifibre()

import casadi.*;

% Add function path
addpath('..\CasadiFunctions\');

% Define parameters
FMo_in = 1;
lMo_in = 1;
lTs_in = 1;
alphao_in = 0;
vMmax_in = [5 10];
load('Fvparam.mat', 'Fvparam');
load('Fpparam.mat', 'Fpparam');
load('Faparam.mat', 'Faparam');
aTendon = 35;
shift = 1;
MuscMoAsmp = false;
d = 0;
stiffness_shift = 1;
stiffness_scale = 1;
strength = 1;
per_fibre = [0.5 0.5];
n_fibre = 2;
N_muscles = 1;

% Function for Hill-equilibrium
FTtilde     = SX.sym('FTtilde', N_muscles); % Normalized tendon forces
a           = SX.sym('a', N_muscles, n_fibre); % Muscle activations
dFTtilde    = SX.sym('dFTtilde', N_muscles); % Time derivative tendon forces
lMT         = SX.sym('lMT', N_muscles); % Muscle-tendon lengths
vMT         = SX.sym('vMT', N_muscles); % Muscle-tendon velocities
tension_SX  = SX.sym('tension', N_muscles); % Tensions
% atendon_SX  = SX.sym('atendon',NMuscle); % Tendon stiffness
% shift_SX    = SX.sym('shift',NMuscle); % shift curve
Hilldiff    = SX(N_muscles, 1); % Hill-equilibrium
FT          = SX(N_muscles, 1); % Tendon forces
Fce         = SX(N_muscles, 1); % Contractile element forces
Fiso        = SX(N_muscles, 1); % Normalized forces from force-length curve
vMmax       = SX(N_muscles, n_fibre); % Maximum contraction velocities
massM       = SX(N_muscles, 1); % Muscle mass
Fpass       = SX(N_muscles, 1); % Passive element forces

% Call function
[err, FT, Fce, Fpass, Fiso, vMmax, massM] = ...
    ForceEquilibrium_FtildeState_all_tendon_multifibre(a, FTtilde, dFTtilde, lMT, ...
    vMT, FMo_in, lMo_in, lTs_in,alphao_in, vMmax_in, Fvparam, Fpparam, ...
    Faparam, tension_SX, aTendon, shift, MuscMoAsmp, d, stiffness_shift, ...
    stiffness_scale, strength, per_fibre, n_fibre);

% Define Casadi function
f_ForceEquilibrium_FtildeState_all_tendon_multifibre = ...
    Function('ForceEquilibrium_FtildeState_all_tendon_multifibre', {a, FTtilde, ...
    dFTtilde, lMT, vMT, tension_SX}, {Hilldiff, FT, Fce, Fpass, Fiso, vMmax, massM}, ...
    {'a', 'FTtilde', 'dFTtilde', 'lMT', 'vMT', 'tension_SX'}, ...
    {'Hilldiff', 'FT', 'Fce', 'Fpass', 'Fiso', 'vMmax', 'massM'});

% Define experiment and integrate over time


end
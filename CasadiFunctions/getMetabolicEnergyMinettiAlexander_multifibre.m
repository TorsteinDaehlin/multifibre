function Edot = getMetabolicEnergyMinettiAlexander_multifibre(a, vM, FMo, vMmax, per_fibre, smeta, n_fibre)
% INPUTS:
%          a: activation           (NMuscle, NFibre)
%         vM: shortening speed     (NMusle, 1)        [m/s]
%     FMo_in: max isometric force  (1, 1)             [N]
%   vMmax_in: max shortening speed (1, NFibre)        [m/s]
%  pre_fibre: fibre fraction       (1, NFibre)
%      smeta: metabolic efficiency (1, NFibre)
%    n_fibre: number of fibres     (1, 1)

% OUTPUTS:
%       Edot: metabolic rate (NMuslce, 1)

% Parameters by fibre
for i = 1:n_fibre
    smeta_sep(:, i) = ones(size(a, 1), 1) * smeta(i);
end

% Minetti & Alexander (1997) model parameters
c1 = 0.054;
c2 = 0.506;
c3 = 2.46;
c4 = 1.13;
c5 = 12.8;
c6 = 1.64;

% Normalized metabolic rate is summed over each fibre before being
% multiplied by maximal force to denormalize
Edot = 0;
for i = 1:n_fibre
    vMtilde = vM ./ vMmax(:, i);

    phi = (c1 + c2 .* (vMtilde) + c3 .* (vMtilde.^2))./ ... 
      (1 - c4 .* (vMtilde) + c5 .* ((vMtilde).^2) - c6.*((vMtilde).^3)); 

    Edot = Edot + a(:, i) .* per_fibre(:, i) .* vMmax(:, i) .* phi .* smeta(:, i);
end
Edot = Edot .* FMo;
end

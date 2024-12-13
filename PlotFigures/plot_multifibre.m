function plot_multifibre(result_path)

% Generate colour palette
cmap = {'#000000','#9c0000', '#1631b8'};

% Load ith results
load(result_path,'R','model_info');

% Get name
[~, name, ~] = fileparts(result_path);

% Plot Qs
figure('Name', [name ' - Qs'], 'WindowState', 'maximized');

nQs = length(model_info.ExtFunIO.coord_names.all);

for i = 1:nQs
    subplot(ceil(sqrt(nQs)), ceil(sqrt(nQs)), i);
    
    plot(R.kinematics.Qs(:, i), 'Color', cmap{1}, 'LineWidth', 1);
    title(replace(R.colheaders.coordinates{i}, '_', ' '));
    xlabel('Time (%)');
    ylabel('Generalized coordinate');
end

% Plot ground reaction forces
figure('Name', [name ' - GRF'], 'WindowState', 'maximized');

dims = size(R.ground_reaction.GRF_r,2);
ax_names = {'X', 'Y', 'Z'};
for i = 1:dims
    subplot(3,3,i);
    plot(R.ground_reaction.GRF_r(:,i), 'Color', cmap{1}, 'LineStyle', '-', 'LineWidth', 1);
    hold on;
    plot(R.ground_reaction.GRF_l(:,i), 'Color', cmap{1}, 'LineStyle', '--', 'LineWidth', 1);
    title(['GRF - ' ax_names{i}]);
    xlabel('Time (%)');
    ylabel('Force (N)');
end

for i = 1:dims
    subplot(3,3,3+i);
    plot(R.ground_reaction.GRM_r(:,i), 'Color', cmap{1}, 'LineStyle', '-', 'LineWidth', 1);
    hold on;
    plot(R.ground_reaction.GRM_l(:,i), 'Color', cmap{1}, 'LineStyle', '--', 'LineWidth', 1);
    title(['GRM - ' ax_names{i}]);
    xlabel('Time (%)');
    ylabel('Moment (Nm)');
end

for i = 1:dims
    subplot(3,3,6+i);
    plot(R.ground_reaction.COP_r(:,i), 'Color', cmap{1}, 'LineStyle', '-', 'LineWidth', 1);
    hold on;
    plot(R.ground_reaction.COP_l(:,i), 'Color', cmap{1}, 'LineStyle', '--', 'LineWidth', 1);
    title(['COP - ' ax_names{i}]);
    xlabel('Time (%)');
    ylabel('Position (m)');
end
legend({'Right', 'Left'}, 'Orientation', 'horizontal', 'Position', [0.5,0.02,0.09,0.02]);

% Plot muscle activations
figure('Name', [name ' - Muscle activations'], 'WindowState', 'maximized');

idx = find(contains(R.colheaders.muscles, {'_r'}));
NMuscle = length(idx);

act_idx = reshape(1:size(R.muscles.a, 2), 2, [])';
act_idx = reshape(act_idx(idx,:)', [], 1);
for i = 1:2:(NMuscle*2)-1
    subplot(ceil(sqrt(NMuscle)), ceil(sqrt(NMuscle)), (i+1)/2);
    
    yline(R.S.bounds.activation_all_muscles.lower, 'k:');
    hold on;
    yline(R.S.bounds.activation_all_muscles.upper, 'k:');
    p1(1) = plot(R.muscles.a(:,act_idx(i)),'Color', cmap{2}, 'LineStyle', '-', 'LineWidth', 1);
    p1(2) = plot(R.muscles.a(:,act_idx(i+1)), 'Color', cmap{3}, 'LineStyle', '-', 'LineWidth', 1);
    title(replace(R.colheaders.muscles{act_idx((i+1)/2)}, '_', ' '));
    ylim([0 1]);
    xlabel('Time (%)');
    ylabel('Activation');
end
legend(p1, {'Slow', 'Fast'}, 'Orientation', 'horizontal', 'Position', [0.5,0.02,0.09,0.02])

% Plot fibre velocities
figure('Name', [name ' - Muscle velocity'], 'WindowState', 'maximized');

for i = 1:2:(NMuscle*2)-1
    subplot(ceil(sqrt(NMuscle)), ceil(sqrt(NMuscle)), (i+1)/2);
    
    p3(1) = plot(R.muscles.vMtilde(:,act_idx(i)),'Color', cmap{2}, 'LineStyle', '-', 'LineWidth', 1);
    hold on;
    p3(2) = plot(R.muscles.vMtilde(:,act_idx(i+1)), 'Color', cmap{3}, 'LineStyle', '-', 'LineWidth', 1);
    title(replace(R.colheaders.muscles{act_idx((i+1)/2)}, '_', ' '));
    xlabel('Time (%)');
    ylabel('Normalized fibre velocity');
    ylim([-1 1]);
end
legend(p3, {'Slow', 'Fast'}, 'Orientation', 'horizontal', 'Position', [0.5,0.02,0.09,0.02])

% Plot muscle forces
figure('Name', [name ' - Muscle forces'], 'WindowState', 'maximized');

for i = 1:NMuscle
    subplot(ceil(sqrt(NMuscle)), ceil(sqrt(NMuscle)), i);

    plot(R.muscles.FTtilde(:,idx(i)), 'Color', cmap{1}, 'LineStyle', '-', 'LineWidth', 1);
    title(replace(R.colheaders.muscles{idx(i)}, '_', ' '));
    xlabel('Time (%)');
    ylabel('Normalized force');
    ylim([0 1]);
end

% Plot metabolic cost
figure('Name', [name ' - Metablic cost'], 'WindowState', 'maximized');

for i = 1:NMuscle
    subplot(ceil(sqrt(NMuscle)), ceil(sqrt(NMuscle)), i);

    plot(R.metabolics.Bhargava2004.Edot_gait(:,idx(i)), 'Color', cmap{1}, 'LineStyle', '-', 'LineWidth', 1);
    title(replace(R.colheaders.muscles{idx(i)}, '_', ' '));
    xlabel('Time (%)');
    ylabel('Energy cost');
    ylim([0 20]);
end

end
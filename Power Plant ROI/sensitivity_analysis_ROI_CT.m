% Sensitivity Analysis

h_orig = 100;  % Sample height
params = [3000, 100000, 75, 20, -0.006, 15, 120, 0.12 * 24 * 365, 10, 150, 500, 10000];
param_names = ["k_material", "k_construction", "T_water", "T_ambient", "gradient_air", "a", "b", "cost_per_MW", "v", "c", "Q_min", "cost_operation_meter_per_year"];

[cost_orig, total_investment_orig] = roi_cooling_tower(h_orig, params);

sensitivity_results = zeros(2, length(params));
for i = 1:length(params)
    % +25%
    params_perturb = params;
    params_perturb(i) = params(i) * 1.25;
    [cost_perturb, ~] = roi_cooling_tower(h_orig, params_perturb);
    sensitivity_results(1, i) = (cost_perturb - cost_orig) / cost_orig * 100;
    
    % -25%
    params_perturb = params;
    params_perturb(i) = params(i) * 0.75;
    [cost_perturb, ~] = roi_cooling_tower(h_orig, params_perturb);
    sensitivity_results(2, i) = (cost_perturb - cost_orig) / cost_orig * 100;
end

% Indices of parameters to plot
plot_indices = [1, 2, 8, 9, 12];  % Indices for "k_material", "k_construction", "cost_per_MW", "v", "cost_operation_meter_per_year" respectively

figure;
hold on;
for i = plot_indices
    plot([1, 2], sensitivity_results(:, i), 'o-', 'DisplayName', param_names(i));
end
hold off;
title('Sensitivity of Cost');
xticks([1, 2]);
xticklabels({'+25%', '-25%'});
ylabel('Percentage Change in Cost (%)');
legend('Location', 'eastoutside');

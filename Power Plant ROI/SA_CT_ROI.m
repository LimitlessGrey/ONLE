%%
close
clear
clc
%%
% Define parameters
temperature = 150;  % initial temperature
targetTemperature = 1;  % final temperature
coolingRate = 0.99;  % cooling rate
maxIterations = 5000;  % maximum number of iterations

% Initialize solution
currentSolution = 150;  % start from the middle of the range
[currentCost, total_investment] = cooling_tower_roi_time(currentSolution);

bestSolution = currentSolution;
bestCost = currentCost;

plot_bestCost = zeros(maxIterations, 1);
plot_currentCost = zeros(maxIterations, 1);

% Main loop
iteration = 1;
while temperature > targetTemperature
    for i = 1:maxIterations
        % Create a new solution
        newSolution = generateNeighbor(currentSolution, temperature);
        [newCost, ~] = cooling_tower_roi_time(newSolution);

        % If new solution is better, accept it
        if newCost < currentCost
            currentSolution = newSolution;
            currentCost = newCost;

            % If new solution is the best so far, update best solution
            if currentCost < bestCost
                bestSolution = currentSolution;
                bestCost = currentCost;
                
            end
        % If new solution is not better, maybe accept it
        else
            accept = rand < exp((currentCost - newCost) / temperature)^2;
            if accept
                currentSolution = newSolution;
                currentCost = newCost;
            end
        end
        % Record costs for plotting
        plot_bestCost(iteration) = bestCost;
        plot_currentCost(iteration) = currentCost;
        iteration = iteration + 1;
    end
    % Decrease the temperature
    temperature = coolingRate * temperature;
end

% Display results
fprintf('Best solution: h = %.2f m\n', bestSolution);
fprintf('Total investment: $%.2f\n', total_investment);
fprintf('ROI time: %.2f years\n', bestCost);

% Plot progress
figure;
hold on;
plot(plot_currentCost, 'b');
plot(plot_bestCost, 'r');
xlabel('Iteration');
ylabel('Cost');
legend('Current Cost', 'Best Cost');
hold off;


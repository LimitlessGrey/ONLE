% Define parameters
temperature = 30;  % initial temperature
targetTemperature = 0.01;  % final temperature
coolingRate = 0.95;  % cooling rate
maxIterations = 100;  % maximum number of iterations
 
% Initialize solution
currentSolution = ones(1, 10) * (40 * 10^-3)^2/4 * pi;

currentCost = Function_Area_SA(currentSolution);  % cost of initial solution
 
% Initialize "best" solution
bestSolution = currentSolution;
bestCost = currentCost;

plot_bestCost = zeros(maxIterations, 1);

currentCostValues = [];
bestCostValues = [];

% Main loop
while temperature > targetTemperature
    for i = 1:maxIterations
        % Create a new solution
        newSolution = generateNeighbor(currentSolution);  % random perturbation of the solution
        newCost = Function_Area_SA(newSolution);  % cost of new solution
 
        % If new solution is better, accept it
        if newCost < currentCost
            currentSolution = newSolution;
            currentCost = newCost;
 
            % If new solution is the best found so far, remember it
            if newCost < bestCost
                bestSolution = newSolution;
                bestCost = newCost;
            end
        % If new solution is worse, maybe accept it
        else
            % Compute acceptance probability
            p_accept = exp(-(newCost - currentCost) / temperature);
 
            % Accept worse solution with a probability
            if rand() < p_accept
                currentSolution = newSolution;
                currentCost = newCost;
            end
        end
        currentCostValues = [currentCostValues, currentCost];
        bestCostValues = [bestCostValues, bestCost];

        % Store the Best Cost Value
        plot_bestCost(i) = bestCost;
        % Display Iteration Information
        disp(['Iteration ' num2str(i) ': Best Cost = ' num2str(plot_bestCost(i))]);

        % Display Temperature Information
        disp(['Temperature ' num2str(temperature)])
    end
 
    % Decrease temperature
    temperature = temperature * coolingRate;
end
 
% Display the best solution
disp(['Best solution: x = ', num2str(bestSolution), ', cost = ', num2str(bestCost)]);

%% Results

figure;
plot(currentCostValues, 'b', 'DisplayName', 'Current Cost');
hold on;
plot(bestCostValues, 'r', 'DisplayName', 'Best Cost');
xlabel('Iteration');
ylabel('Objective Function Value');
title('Objective Function Value over Iterations');
legend('Location', 'best');
grid on;

%% Generate a neighboring solution
function newSolution = generateNeighbor(currentSolution)
    % Generate a neighboring solution based on the current solution
    
    % Determine the minimum value allowed for each component
    minValue = ((0.005^2) / (4 * pi));
    
    % Initialize the new solution as a copy of the current solution
    newSolution = currentSolution;
    
    % Iterate over each component and make adjustments
    for i = 1:10
        % Generate a random value within a certain range
        adjustment = randn() * 0.00001; % Adjust the range as needed
        
        % Apply the adjustment to the current component
        newSolution(i) = currentSolution(i) + adjustment;
        
        % Apply the constraint: if the adjusted value is below the minimum, set it to the minimum
        if newSolution(i) < minValue
            newSolution(i) = minValue;
        end
    end
end

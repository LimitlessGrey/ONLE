% Define parameters
temperature = 50;  % initial temperature
targetTemperature = 0.1;  % final temperature
coolingRate = 0.97;  % cooling rate
maxIterations = 500;  % maximum number of iterations
 
% Initialize solution
currentSolution = [0.5,0.5,0.5,0.5];

currentCost = function_evaluate_time(currentSolution);  % cost of initial solution
 
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
        newCost = function_evaluate_time(newSolution);  % cost of new solution
 
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
    % Initialize the new solution as a copy of the current solution
    newSolution = currentSolution;
    
    % Determine the range for perturbation (you can adjust these values)
    perturbationRange = 0.01;
    
    % Apply random perturbation to each component of the vector
    for i = 1:length(currentSolution)
        perturbation = perturbationRange * (rand() - 0.5);  % generate a random perturbation
        newSolution(i) = newSolution(i) + perturbation;  % apply the perturbation
        
        % Ensure the component is greater than the threshold value
        if newSolution(i) <= 0.246609125
            newSolution(i) = 0.246609125 + abs(perturbation);  % adjust the value
        end
    end
end
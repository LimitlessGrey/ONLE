function newSolution = generateNeighbor(currentSolution, temperature)
% Generate a neighboring solution based on the current solution

% Determine the minimum and maximum value for h
minValue = 100;
maxValue = 200;

% Generate a random value within a range proportional to the current temperature
adjustment = randn() * (10 * temperature / 1500); % Adjust the range as needed
% Limit the adjustment so that it can't exceed the range
adjustment = max(min(adjustment, maxValue - currentSolution), minValue - currentSolution);

% Apply the adjustment to the current solution
newSolution = currentSolution + adjustment;

% Apply the constraint: if the adjusted value is outside the bounds, clip it
if newSolution < minValue
    newSolution = minValue;
elseif newSolution > maxValue
    newSolution = maxValue;
end
end

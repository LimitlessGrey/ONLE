function time_total = function_evaluate_time(x)
%% Model creation

data = readmatrix('power_curve.txt');
data(:, 2) = data(:, 2) * 735.5; % Convert [HP] to [W]
data(:, 1) = data(:, 1) * 2*pi/60; % Convert [RPM] to [rad/s]

eff = 0.9; % Transmission efficiency
diff = 0.345265; % Differential ratio
R = 0.3268; % Tyre radius [m]

% Resisting power parameters
A = 50;
B = 0.5272;

x = [0.246609125, x(1), x(2), x(3), x(4), 1.14678899]; % Gear Ratios

num_points = size(data, 1); % Get number of rows of data matrix
num_gears = size(x, 2); % Get number of gear ratios

v = zeros(num_points, num_gears); % Initialize velocity matrix
Pe = zeros(num_points, num_gears); % Initialize power matrix

for j=1:num_gears
    for i=1:num_points
        v(i, j) = R * diff * x(j) * data(i, 1); % [m/s]
        Pn = A * v(i, j) + B * v(i, j)^3; % Resisting power [W]
        Pw = data(i, 2) * eff; % Power to the wheels [W]

        Pe(i, j) = Pw - Pn; % Excess power [W]
    end
end

%% Time evaluation

m = [2640, 2475, 2310, 2145, 1980, 1815]; % Apparent mass vector
t_gc = 0.5; % Gear change time
target_velocity = 300/3.6;

[~, max_power_index] = max(data(:, 2)); % Store max power and its index

time_integral = zeros(1, num_gears);

j = 1;
while j <= num_gears
    v_gear = v(1:max_power_index, j);
    Pe_gear = Pe(1:max_power_index, j);
    if j == 1
        % Calculate the time integral using the trapezoidal rule
        time_integral(j) = trapz(v_gear, m(j) .* v_gear ./ Pe_gear);
    else
        % Get the previous gear's velocity
        v_previous_gear = v(1:max_power_index, j-1);
        
        % Find the index of the closest velocity in the next gear's velocity vector
        [~, index] = min(abs(v_gear - v_previous_gear(max_power_index)));

        % Check if the stopping velocity is reached in the current gear
        if v_gear(max_power_index) >= target_velocity
            % Find the index where the stopping velocity is reached
            stopping_index = find(v_gear >= target_velocity, 1);
            % Calculate the time integral up to the stopping velocity
            time_integral(j) = trapz(v_gear(index-1:stopping_index), m(j) .* v_gear(index-1:stopping_index) ./ Pe_gear(index-1:stopping_index));
            break;  % Exit the loop since the stopping criterion is met
        elseif v_previous_gear(max_power_index) > max(v_gear)
            % Skip this gear since the velocity is not found
            j = j + 1;
            disp('skipped');
            g = -(x(j) - x(j-1));
            r = 10;
            time_integral(j) = time_integral(j) + r * max(0, g)^2;
            continue;
        else
            % Calculate the time integral for the current gear using all available velocities
            time_integral(j) = trapz(v_gear(index-1:max_power_index), m(j) .* v_gear(index-1:max_power_index) ./ Pe_gear(index-1:max_power_index));
        end
    end
    j = j + 1;
end
time_total = sum(time_integral) + (j-1)*t_gc;
end
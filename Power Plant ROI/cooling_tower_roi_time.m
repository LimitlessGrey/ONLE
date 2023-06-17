%% Objective Function
function [cost, total_investment] = cooling_tower_roi_time(h)
% Parameters
k_material = 3000;  % $/m^2
k_construction = 100000;  % $/m
T_water = 75;  % C
T_ambient = 20;  % C
gradient_air = -0.006;  % C/m
a = 15;  % W/(m^2K)
b = 120;  % W/(m^2K)
cost_per_MW = 0.12 * 24 * 365;  % $/MW/year, adjusted price
v = 10;  % m/s
c = 150;  % constant for energy production calculation, adjust as needed
Q_min = 500;  % MW, minimum heat transfer requirement

% Operation cost per meter per year
cost_operation_meter_per_year = 10000;  % $/m/year

% Geometry
D1 = 60 + 0.4 * h;  % diameter at the bottom
D2 = 40 + 0.2 * h;  % diameter at the top

% Area
A = pi * (D1 + D2) * sqrt((D1 - D2)^2 / 4 + h^2);

% Cost of materials and construction
C_material = k_material * A;
C_construction = k_construction * h;

% Temperature at height h
T_air_h = T_ambient + h^2 * gradient_air;

% Heat transfer coefficient
hc = a * v + b;

% Actual heat transfer
Q_actual = hc * A * (T_water - T_air_h) / 1e6;  % MW

% Energy production
energy_production = c * h;  % Adjust the function as needed

% ROI Calculation
revenue_per_year = energy_production * cost_per_MW;  % $/year
total_investment = C_material + C_construction;  % $

% Cost of operation (annual)
if Q_actual < Q_min
   Cost_of_Operation_per_year = h * cost_operation_meter_per_year + 5e10;  
else 
   Cost_of_Operation_per_year = h * cost_operation_meter_per_year;
end

% Total annual Profit
total_annual_profit = revenue_per_year - Cost_of_Operation_per_year;

% Time for ROI
ROI_time = total_investment / total_annual_profit;  % years

% Total cost is time for ROI
cost = ROI_time;

% Debugging print statements
% fprintf('Height: %.2f m\n', h);
% fprintf('Q_actual: %.2f MW\n', Q_actual);
% fprintf('Energy production: %.2f MW\n', energy_production);
% fprintf('Revenue per year: $%.2f\n', revenue_per_year);
% fprintf('Total investment: $%.2f\n', total_investment);
% fprintf('Cost of Operation per year: $%.2f\n', Cost_of_Operation_per_year);
% fprintf('Total annual profit: $%.2f\n', total_annual_profit);
% fprintf('ROI time: %.2f years\n', ROI_time);
end

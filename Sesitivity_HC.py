#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Define the objective function (adapt to your specific problem)
def objective_function(a, *params):
    # Calculate the fuel consumption based on the acceleration
    # fuel consumption (a)
    return fuel_consumption

# Define the constraints (adapt to your specific problem)
constraints = (
    # Velocity constraint (Max Vel = 120kmh)
    Vmax = 120,
    # Power(force) constraint (Still trying to understand how to correctly define power from accelration)
    Ft = 
    # Time constraint (Ma)
    MaxT = 3
)

# Define the acceleration values (might be Binary)
acceleration_bounds = [(a_min, a_max) for _ in range(N)]

# Define a range of acceleration values for the sensitivity analysis
a_test_values = np.linspace(a_min, a_max, num_test_values)

# Store the optimal fuel consumption results
optimal_fuel_consumption_results = []

# Perform sensitivity analysis for each acceleration value
for a_test in a_test_values:
    # Update the constraints with the test acceleration value
    # ...

    # Solve the optimization problem
    result = minimize(
        objective_function,
        x0=np.zeros(N),
        bounds=acceleration_bounds,
        constraints=constraints,
        args=(params)
    )

    # Append if the optimization went OK
    if result.success:
        optimal_fuel_consumption_results.append(result.fun)
    else:
        optimal_fuel_consumption_results.append(None)

# Analyze the results Plot the 

plt.plot(a_test_values, optimal_fuel_consumption_results)
plt.xlabel('Acceleration')
plt.ylabel('Optimal Fuel Consumption')
plt.title('Sensitivity Analysis of Acceleration on Fuel Consumption')
plt.show()
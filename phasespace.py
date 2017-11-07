"""
performs the sympletic Euler method for a mass on a spring
uses previously defined methods from euler.py (explicit, implicit, analytic)
calculate normalized energy for all methods using function definied in euler.py
plots phase-space geometries from sympletic and other methods
plot energies of sympletic and analytic solutions
"""

import numpy as np
import matplotlib.pyplot as plt
import euler

# calculates mass on spring using sympletic Euler method
def sympletic(total_steps, step_size, x_0, v_0):
    timesteps = np.arange(0, total_steps, step_size)
    x_values = np.zeros(len(timesteps))
    v_values = np.zeros(len(timesteps))
    x_values[0] = x_0
    v_values[0] = v_0
    for i in range(len(timesteps)-1):
        x_values[i+1] = x_values[i] + (step_size * v_values[i])
        v_values[i+1] = v_values[i] - (step_size * x_values[i+1])
    return x_values, v_values

# set initial parameters
total_steps = 10
step_size = .1
x_0 = 1
v_0 = 0

# create array of t values
timesteps = np.arange(0, total_steps, step_size)

# use previously defined explicit and implicit Euler methods to calculate x and v values
x_explicit, v_explicit = euler.explicit(total_steps, step_size, x_0, v_0)
x_implicit, v_implicit = euler.implicit(total_steps, step_size, x_0, v_0)

# calculate analytic solution (previously defined)
x_analytic, v_analytic = euler.analytic(total_steps, step_size, x_0, v_0)

# calculate sympletic values for x and v
x_sympletic, v_sympletic = sympletic(total_steps, step_size, x_0, v_0)

# calculate energies for sympletic and analytic solutions
energy_analytic = euler.energy(total_steps, step_size, x_analytic, v_analytic)
energy_sympletic = euler.energy(total_steps, step_size, x_sympletic, v_sympletic)

# plot phase-space geometries; uncomment desired lines
#plt.plot(x_explicit, v_explicit, '-', label='explicit')
#plt.plot(x_implicit, v_implicit, '-', label='implicit')
#plt.plot(x_analytic, v_analytic, '-', label='analytic')
#plt.plot(x_sympletic, v_sympletic, '-', label='sympletic')

# plot energies of sympletic and analytic solutions (make sure all above plots are commented out!)
#plt.plot(timesteps, energy_analytic, '-', label='analytic')
#plt.plot(timesteps, energy_sympletic, '-', label='sympletic')

plt.legend()
plt.show()
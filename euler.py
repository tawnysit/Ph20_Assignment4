"""
this program performs the explicit or implicit Euler integration for a mass on a spring
it also can calculate the analytic solution to the problem
using these solutions, global and truncation error and normalized energy can be calculated
can choose which values to plot by uncommenting lines at the end
"""

import numpy as np
import matplotlib.pyplot as plt

# calculates mass on spring using explicit Euler method
def explicit(total_steps, step_size, x_0, v_0):
    timesteps = np.arange(0, total_steps, step_size)
    x_values = np.zeros(len(timesteps))
    v_values = np.zeros(len(timesteps))
    x_values[0] = x_0
    v_values[0] = v_0
    for i in range(len(timesteps)-1):
        x_values[i+1] = x_values[i] + (step_size * v_values[i])
        v_values[i+1] = v_values[i] - (step_size * x_values[i])
    return x_values, v_values

# calculates mass on spring using implicit Euler method
def implicit(total_steps, step_size, x_0, v_0):
    timesteps = np.arange(0, total_steps, step_size)
    x_values = np.zeros(len(timesteps))
    v_values = np.zeros(len(timesteps))
    x_values[0] = x_0
    v_values[0] = v_0
    for i in range (len(timesteps)-1):
        x_values[i+1] = (x_values[i] + (step_size * v_values[i])) / (step_size**2 + 1)
        v_values[i+1] = (v_values[i] - (step_size * x_values[i])) / (step_size**2 + 1)
    return x_values, v_values

# analytical solution to mass on a spring
def analytic(total_steps, step_size, x_0, v_0):
    timesteps = np.arange(0, total_steps, step_size)
    x_values = np.zeros(len(timesteps))
    v_values = np.zeros(len(timesteps))
    x_values[0] = x_0
    v_values[0] = v_0
    for i in range(len(timesteps)):
        x_values[i] = (x_0 * np.cos(timesteps[i])) + (v_0 * np.sin(timesteps[i]))
        v_values[i] = (-1 * x_0 * np.sin(timesteps[i])) + (v_0 * np.cos(timesteps[i]))
    return x_values, v_values

# computes global error based on analytical solution and numerical integration
def global_error(x_analytic, v_analytic, x, v):   
    error_x = np.zeros(len(timesteps))
    error_v = np.zeros(len(timesteps))
    for i in range(len(timesteps)):
        error_x[i] = x_analytic[i] - x[i]
        error_v[i] = v_analytic[i] - v[i]
    return error_x, error_v

# computes truncation error for a set of step sizes for the explicit Euler method
def truncation_error(h_0, total_time, x_0, v_0):
    step_sizes = np.zeros(5)
    step_sizes[0] = h_0
    for i in [1,2,3,4]:
        step_sizes[i] = h_0 / (2**i)
    t_error = np.zeros(len(step_sizes))
    j = 0
    for h in step_sizes:
        x, v = explicit(total_time, h, x_0, v_0)
        x_analytic, v_analytic = analytic(total_time, h, x_0, v_0)
        error_x = global_error(x_analytic, v_analytic, x, v)[0]
        t_error[j] = np.amax(np.absolute(error_x))
        j += 1
    return step_sizes, t_error

# computes normalized energy as a function of time
def energy(total_steps, step_size, x, v):
    timesteps = np.arange(0, total_steps, step_size)
    energy = np.zeros(len(timesteps))
    for i in range(len(timesteps)):
        energy[i] = x[i]**2 + v[i]**2
    return energy

#set initial parameters
total_steps = 20
step_size = .1
x_0 = 1
v_0 = 0

# creates array of t values, same as in functions defined above
timesteps = np.arange(0, total_steps, step_size)

# calculate explicit and implicit values for x and v
x_explicit, v_explicit = explicit(total_steps, step_size, x_0, v_0)
x_implicit, v_implicit = implicit(total_steps, step_size, x_0, v_0)

# calculate analytic solution
x_analytic, v_analytic = analytic(total_steps, step_size, x_0, v_0)

# calculate global error and energy for both explicit and implicit integration
error_x_explicit, error_v_explicit = global_error(x_analytic, v_analytic, x_explicit, v_explicit)
error_x_implicit, error_v_implicit = global_error(x_analytic, v_analytic, x_implicit, v_implicit)
energy_explicit = energy(total_steps, step_size, x_explicit, v_explicit)
energy_implicit = energy(total_steps, step_size, x_implicit, v_implicit)

#calculate truncation error
h_values, truncation_error = truncation_error(step_size, total_steps, x_0, v_0)

# uncomment line(s) that you wish to plot:

#plt.plot(timesteps, x_explicit, '-', label='x (explicit)') # plot numerical x values from explicit integration
#plt.plot(timesteps, v_explicit, '-', label='v (explicit)') # plot numerical v values from explicit integration
#plt.plot(timesteps, x_implicit, '-', label='x (implicit)') # plot numerical x values from implicit integration
#plt.plot(timesteps, v_implicit, '-', label='v (implicit)') # plot numerical v values from implicit integration
#plt.plot(timesteps, x_analytic, '-', label='x analytic') # plot analytic x values
#plt.plot(timesteps, v_analytic, '-', label='v analytic') # plot analytic v values
#plt.plot(timesteps, error_x_explicit, '-', label='explicit x error') # plot global error in x from explicit integration
#plt.plot(timesteps, error_v_explicit, '-', label='explicit v error') # plot global error in v from explicit integration
#plt.plot(timesteps, error_x_implicit, '-', label='implicit x error') # plot global error in x from implicit integration
#plt.plot(timesteps, error_v_implicit, '-', label='implicit v error') # plot global error in v from implicit integration
#plt.plot(h_values, truncation_error, 'o') # plot truncation error as function of step size
#plt.plot(timesteps, energy_explicit, '-') # plot energy of numerical integration from explicit integration
#plt.plot(timesteps, energy_implicit, '-') # plot energy of numerical integration from implicit integration

plt.legend()
plt.show()
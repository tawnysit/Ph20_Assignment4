import numpy as np
import matplotlib.pyplot as plt
import sys

# calculates mass on spring using explicit Euler method
def explicit(start, end, step_size, x_0, v_0):
    steps = (end-start)/step_size
    timesteps = np.linspace(start, end, steps+1)
    x_values = np.zeros(len(timesteps))
    v_values = np.zeros(len(timesteps))
    x_values[0] = x_0
    v_values[0] = v_0
    for i in range(len(timesteps)-1):
        x_values[i+1] = x_values[i] + (step_size * v_values[i])
        v_values[i+1] = v_values[i] - (step_size * x_values[i])
    return x_values, v_values

# calculates mass on spring using implicit Euler method
def implicit(start, end, step_size, x_0, v_0):
    steps = (end-start)/step_size
    timesteps = np.linspace(start, end, steps+1)
    x_values = np.zeros(len(timesteps))
    v_values = np.zeros(len(timesteps))
    x_values[0] = x_0
    v_values[0] = v_0
    for i in range (len(timesteps)-1):
        x_values[i+1] = (x_values[i] + (step_size * v_values[i])) / (step_size**2 + 1)
        v_values[i+1] = (v_values[i] - (step_size * x_values[i])) / (step_size**2 + 1)
    return x_values, v_values

# analytical solution to mass on a spring
def analytic(start, end, step_size, x_0, v_0):
    steps = (end-start)/step_size
    timesteps = np.linspace(start, end, steps+1)
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
        x, v = explicit(start, total_time, h, x_0, v_0)
        x_analytic, v_analytic = analytic(start, total_time, h, x_0, v_0)
        error_x = global_error(x_analytic, v_analytic, x, v)[0]
        t_error[j] = np.amax(np.absolute(error_x))
        j += 1
    return step_sizes, t_error

# computes normalized energy as a function of time
def energy(start, end, step_size, x, v):
    steps = (end-start)/step_size
    timesteps = np.linspace(start, end, steps+1)
    energy = np.zeros(len(timesteps))
    for i in range(len(timesteps)):
        energy[i] = x[i]**2 + v[i]**2
    return energy

# calculates mass on spring using sympletic Euler method
def sympletic(start, end, step_size, x_0, v_0):
    steps = (end-start)/step_size
    timesteps = np.linspace(start, end, steps+1)
    x_values = np.zeros(len(timesteps))
    v_values = np.zeros(len(timesteps))
    x_values[0] = x_0
    v_values[0] = v_0
    for i in range(len(timesteps)-1):
        x_values[i+1] = x_values[i] + (step_size * v_values[i])
        v_values[i+1] = v_values[i] - (step_size * x_values[i+1])
    return x_values, v_values

#set initial parameters
start = 0
end = 20
step_size = .1
x_0 = 1
v_0 = 0

# creates array of t values, same as in functions defined above
steps = (end-start)/step_size
timesteps = np.linspace(start, end, steps+1)

# calculate explicit, implicit, sympletic values for x and v
x_explicit, v_explicit = explicit(start, end, step_size, x_0, v_0)
x_implicit, v_implicit = implicit(start, end, step_size, x_0, v_0)
x_sympletic, v_sympletic = sympletic(start, end, step_size, x_0, v_0)

# calculate analytic solution
x_analytic, v_analytic = analytic(start, end, step_size, x_0, v_0)

# calculate global error and energy for both explicit and implicit integration
error_x_explicit, error_v_explicit = global_error(x_analytic, v_analytic, x_explicit, v_explicit)
error_x_implicit, error_v_implicit = global_error(x_analytic, v_analytic, x_implicit, v_implicit)
energy_explicit = energy(start, end, step_size, x_explicit, v_explicit)
energy_implicit = energy(start, end, step_size, x_implicit, v_implicit)

# calculate energies for sympletic and analytic solutions
energy_analytic = energy(start, end, step_size, x_analytic, v_analytic)
energy_sympletic = energy(start, end, step_size, x_sympletic, v_sympletic)

#calculate truncation error
h_values, truncation_error = truncation_error(step_size, end, x_0, v_0)

# plot and save x and v values from explicit integration
plt.plot(timesteps, x_explicit, '-', label='x (explicit)')
plt.plot(timesteps, v_explicit, '-', label='v (explicit)')
plt.legend()
plt.savefig("explicit_solution.jpg")
plt.close()

# plot and save x and v values from implicit integration
plt.plot(timesteps, x_implicit, '-', label='x (implicit)')
plt.plot(timesteps, v_implicit, '-', label='v (implicit)')
plt.legend()
plt.savefig("implicit_solution.jpg")
plt.close()

# plot and save x and v values from analytic solution
plt.plot(timesteps, x_analytic, '-', label='x (analytic)')
plt.plot(timesteps, v_analytic, '-', label='v (analytic)')
plt.legend()
plt.savefig("analytic_solution.jpg")
plt.close()

# plot and save explicit vs implicit solution
plt.plot(timesteps, x_explicit, '-', label='x (explicit)')
plt.plot(timesteps, v_explicit, '-', label='v (explicit)')
plt.plot(timesteps, x_implicit, '-', label='x (implicit)')
plt.plot(timesteps, v_implicit, '-', label='v (implicit)')
plt.legend()
plt.savefig("implicit_vs_explicit.jpg")
plt.close()

# plot and save analytic vs explicit solution
plt.plot(timesteps, x_analytic, '-', label='x (analytic)')
plt.plot(timesteps, v_analytic, '-', label='v (analytic)')
plt.plot(timesteps, x_explicit, '-', label='x (explicit)')
plt.plot(timesteps, v_explicit, '-', label='v (explicit)')
plt.legend()
plt.savefig("analytic_vs_explicit.jpg")
plt.close()

# plot and save analytic vs implicit solution
plt.plot(timesteps, x_analytic, '-', label='x (analytic)')
plt.plot(timesteps, v_analytic, '-', label='v (analytic)')
plt.plot(timesteps, x_implicit, '-', label='x (implicit)')
plt.plot(timesteps, v_implicit, '-', label='v (implicit)')
plt.legend()
plt.savefig("analytic_vs_implicit.jpg")
plt.close()

# plot and save errors from explicit integration
plt.plot(timesteps, error_x_explicit, '-', label='explicit x error')
plt.plot(timesteps, error_v_explicit, '-', label='explicit v error')
plt.legend()
plt.savefig("explicit_error.jpg")
plt.close()

# plot and save errors from implicit integration
plt.plot(timesteps, error_x_implicit, '-', label='implicit x error')
plt.plot(timesteps, error_v_implicit, '-', label='implicit v error')
plt.legend()
plt.savefig("implicit_error.jpg")
plt.close()

# plot and save errors from implicit and explicit integration
plt.plot(timesteps, error_x_explicit, '-', label='explicit x error')
plt.plot(timesteps, error_v_explicit, '-', label='explicit v error')
plt.plot(timesteps, error_x_implicit, '-', label='implicit x error')
plt.plot(timesteps, error_v_implicit, '-', label='implicit v error')
plt.legend()
plt.savefig("implicit_vs_explicit_error.jpg")
plt.close()

# plot and save truncation error
plt.plot(h_values, truncation_error, 'o')
plt.savefig("truncation_error.jpg")
plt.close()

# plot energy from explicit integration
plt.plot(timesteps, energy_explicit, '-')
plt.savefig("explicit_energy.jpg")
plt.close()

# plot energy from implicit integration
plt.plot(timesteps, energy_implicit, '-')
plt.savefig("implicit_energy.jpg")
plt.close()

# implicit, explicit, analytic phasespace
plt.plot(x_explicit, v_explicit, '-', label='explicit')
plt.plot(x_implicit, v_implicit, '-', label='implicit')
plt.plot(x_analytic, v_analytic, '-', label='analytic')
plt.legend()
plt.savefig("implicit_vs_explicit_vs_analytic_phasespace.jpg")
plt.close()

# sympletic phasespace
plt.plot(x_sympletic, v_sympletic, '-', label='sympletic')
plt.legend()
plt.savefig("sympletic_phasespace.jpg")
plt.close()

# analytic phasespace
plt.plot(x_analytic, v_analytic, '-', label='analytic')
plt.legend()
plt.savefig("analytic_phasespace.jpg")
plt.close()

# analytic vs sympletic phasespace
plt.plot(x_sympletic, v_sympletic, '-', label='sympletic')
plt.plot(x_analytic, v_analytic, '-', label='analytic')
plt.legend()
plt.savefig("analytic_vs_sympletic_phasespace.jpg")
plt.close()

# plot energies of sympletic and analytic solutions
plt.plot(timesteps, energy_analytic, '-', label='analytic')
plt.plot(timesteps, energy_sympletic, '-', label='sympletic')
plt.legend()
plt.savefig("analytic_vs_sympletic_energy.jpg")
plt.close()
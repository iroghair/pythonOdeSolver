import numpy as np
from scipy.integrate import odeint
from odeintegrator import myodeint
import matplotlib.pyplot as plt


def myODEfunc(y,t,k=0.2):
    dydt = -k*y

    return dydt


# Initial  conditions 
init0 = 1

# Time steps
steps = np.arange(0, 10, 0.1) # Move out of iteration loop

# Other parameters
k = 0.3

# ODE Solver
sol = myodeint(myODEfunc, init0, steps)

plt.plot(steps, sol)
# plt.legend(loc='best')
plt.xlabel('Time [s]')
plt.ylabel('Concentration')
plt.title('ODE Solution')
plt.grid()
plt.show()
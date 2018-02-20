import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def myODEfunc(y,t,k):
    dydt = -k*y

    return dydt


# Initial  conditions 
init0 = 1

# Time steps
steps = np.arange(0, 10, 0.1) # Move out of iteration loop

# Other parameters
k = 0.3

# ODE Solver
sol = odeint(myODEfunc, init0, steps, args=(k,))

plt.plot(steps, sol)
# plt.legend(loc='best')
plt.xlabel('Time [s]')
plt.ylabel('Concentration')
plt.title('ODE Solution')
plt.grid()
plt.show()
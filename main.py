import numpy as np
from scipy.integrate import odeint
from odeintegrator import myodeint
import matplotlib.pyplot as plt


def myODEfunc(y,t,k=0.1):
    dydt = -k*y

    return dydt


# Initial  conditions 
init0 = 1

# Time steps
steps = np.linspace(0, 10, num=11) # Move out of iteration loop

# Other parameters
k = 0.3

# ODE Solver
span,sol = myodeint(myODEfunc, init0, steps, solver='Fehlberg')

plt.plot(span, sol, marker="o")
# plt.legend(loc='best')
plt.xlabel('Time [s]')
plt.ylabel('Concentration')
plt.title('ODE Solution')
plt.grid()
plt.show()
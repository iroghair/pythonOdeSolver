import numpy as np
import collections
import sys

def myodeint(func, init, span, solver='Euler', 
    N=1000, **kwargs):
    
    if len(span) == 2:
        # Set dt if tspan of 2 is given
        print("Span is of length 2, defaulting to ", N, " steps.")
        span = np.linspace(span[0],span[1],num=N+1)

    try:
        dt = span[1]-span[0]
        y = init
    except IndexError:
        print("Error! The time span should contain either the start and end times, or all time steps desired.")
        sys.exit(0)

    sol = np.zeros( [len(span), len([init])] )
    sol[0,:] = init

    if solver == 'Euler':
        for t in range(1,len(span)):
            tdt = t+dt
            k1 = func(y,t)
            y = y+dt*k1
            sol[t,:] = y
    elif solver == 'Midpoint':
        for t in range(1,len(span)):
            dh = dt/2
            tdh = t+dh
            k1 = func(y,t)
            q1 = y+dh*k1
            k2 = func(q1,tdh)
            y = y+dt*k2
            sol[t,:] = y
    elif solver == 'RK2':
        for t in range(1,len(span)):
            dh = dt/2
            tdh = t+dh
            tdt = t+dt
            k1 = func(y,t)
            q1 = y+dt*k1
            k2 = func(q1,tdt)
            y = y+dt*0.5*(k1+k2)
            sol[t,:] = y
    elif solver == 'RK4':
        for t in range(1,len(span)):
            dh = dt/2
            tdh = t+dh
            tdt = t+dt
            k1 = func(y,t)
            q1 = y+0.5*dt*k1
            k2 = func(q1,tdh)
            q2 = y+0.5*dt*k2
            k3 = func(q2,tdh)
            q3 = y+dt*k3
            k4 = func(q3,tdt)
            y = y+dt*(k1+2*k2+2*k3+k4)/6
            sol[t,:] = y
#     case 'Verlet'
#         % Store initial values
#         x = init(1);
#         v = init(2);
#         % First take an Euler step
#         k1 = func(t,y,par);
#         y = y+dt*k1;
#         x(2) = y(1);
#         v(2) = y(2);
#         % Set the first accelerations
#         a(1) = 0;
#         a(2) = k1(2);
#         for t = 3:length(span);
#             x(t) = x(t-1) + v(t-1)*dt + 1/6*(4*a(t-1) - a(t-2))*dt^2;
#             k = func(t,[x(t) v(t-1)],par);
#             a(t) = k(2); % store acceleration again
#             v(t) = v(t-1) + 1/6 * (2*a(t) + 5*a(t-1) - a(t-2))*dt;
#         end
#         sol = [x; v]';
    elif solver == 'Implicit':
        dh = dt/2
        for t in range(1,len(span)):
            # Numerical derivative via Euler time step: d/dy*(dy/dx)
            k1 = func(y,t)
            y2 = y+dt*k1
            k2 = func(y,t+dt)
            dy = np.asarray(y2-y)
            dy[dy==0] = 1e-12
            dfdy = [(k2-k1)] / dy
            # Compute next time step using implicit
            y = y+dt*(1/(1-dh*dfdy))*k1
            sol[t,:] = y
    elif solver == 'Fehlberg':
        sol = np.zeros( [len(span), 2])
        sol[0,:] = init

        c1 = 0
        c2 = 1/5
        c3 = 3/10
        c4 = 4/5
        c5 = 8/9
        c6 = 1
        c7 = 1

        a21 = 1/5
        a31 = 3/40
        a32 = 9/40
        a41 = 44/45
        a42 = -56/15
        a43 = 32/9
        a51 = 19372/6561
        a52 = -25360/2187
        a53 = 64448/6561
        a54 = -212/729
        a61 = 9017/3168
        a62 = -355/33
        a63 = 46732/5247
        a64 = 49/176
        a65 = -5103/18656

        b1 = 35/384
        b2 = 0
        b3 = 500/1113
        b4 = 125/192
        b5 = -2187/6784
        b6 = 11/84

        b1s = 5179/57600
        b2s = 0
        b3s = 7571/16695
        b4s = 393/640
        b5s = -92097/339200
        b6s = 187/2100

        for t in range(1,len(span)):
            k1 = func(y,t)
            k2 = func(y+a21*k1, t+c2*dt)
            k3 = func(y+a31*k1 + a32*k2, t+c3*dt)
            k4 = func(y+a41*k1 + a42*k2 + a43*k3, t+c4*dt)
            k5 = func(y+a51*k1 + a52*k2 + a53*k3 + a54*k4, t+c5*dt)
            k6 = func(y+a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5, t+c6*dt)

            sol[t,:] = [y + dt*(b1*k1  + b2*k2  + b3*k3  + b4*k4  + b5*k5  + b6*k6 ),
                        y + dt*(b1s*k1 + b2s*k2 + b3s*k3 + b4s*k4 + b5s*k5 + b6s*k6)]

            y = sol[t,1]

#     otherwise
#         disp('Method not supported - Use Euler, Midpoint, RK2, RK4, Verlet');
#         span = []; 
#         sol = [];
#         return
# end

# end

    return span, sol

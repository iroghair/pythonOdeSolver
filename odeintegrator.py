import numpy as np

def myodeint(func, init, span, solver='Euler',N=1000):
    
    if len(span) == 2:
        # Set dt if tspan of 2 is given
        print("Span is of length 2, defaulting to ", N, " steps.")
        span = np.range(span[0],span[1],1+N)

    
    dt = span[1]-span[0]

    t = span[0]
    y = init

    sol = np.zeros( [len(span), len([init])] )
    sol[0,:] = init

    if solver == 'Euler':
        for t in range(2,len(span)):
            tdt = t+dt;
            k1 = func(y,t);
            y = y+dt*k1;
            sol[t,:] = y;
#     case 'Midpoint'
#         for t = 2:length(span);
#             dh = dt/2;
#             tdh = t+dh;
#             k1 = func(t,y,par);
#             q1 = y+dh*k1;
#             k2 = func(tdh,q1,par);
#             y = y+dt*k2;
#             sol(t,:) = y;
#         end
#     case 'RK2'
#         for t = 2:length(span);
#             dh = dt/2;
#             tdh = t+dh;
#             tdt = t+dt;
#             k1 = func(t,y,par);
#             q1 = y+dt*k1;
#             k2 = func(tdt,q1,par);
#             y = y+dt*0.5*(k1+k2);
#             sol(t,:) = y;
#         end
#     case 'RK4'
#         for t = 2:length(span);
#             dh = dt/2;
#             tdh = t+dh;
#             tdt = t+dt;
#             k1 = func(t,y,par);
#             q1 = y+0.5*dt*k1;
#             k2 = func(tdh,q1,par);
#             q2 = y+0.5*dt*k2;
#             k3 = func(tdh,q2,par);
#             q3 = y+dt*k3;
#             k4 = func(tdt,q3,par);
#             y = y+dt*(k1+2*k2+2*k3+k4)/6;
#             sol(t,:) = y;
#         end
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
#     case 'Implicit'
#         dh = dt/2;
#         for t = 2:length(span);
#             % Numerical derivative via Euler time step: d/dy*(dy/dx) 
#             k1 = func(t,y,par);
#             y2 = y+dt*k1;
#             k2 = func(t+dt,y2,par);
#             dy = (y2-y);
#             dy(dy==0) = 1e-6;
#             dfdy = (k2-k1)./dy;
#             % Compute next time step using implicit
#             y = y+dt*(1/(1-dh*dfdy))*k1;
#             sol(t,:) = y;
#         end
#     otherwise
#         disp('Method not supported - Use Euler, Midpoint, RK2, RK4, Verlet');
#         span = []; 
#         sol = [];
#         return
# end

# end

    return sol, span

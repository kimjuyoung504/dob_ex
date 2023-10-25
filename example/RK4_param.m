function y = RK4_param(func, t, y, u, dt, param)

k_1 = func(t     , y           , u, param);
k_2 = func(t+dt/2, y+(dt/2)*k_1, u, param);
k_3 = func(t+dt/2, y+(dt/2)*k_2, u, param);
k_4 = func(t+dt  , y+    dt*k_3, u, param);

y = y + (dt/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
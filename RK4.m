function y = RK4(func, t, y, u, dt)

k_1 = func(t     , y           , u);
k_2 = func(t+dt/2, y+(dt/2)*k_1, u);
k_3 = func(t+dt/2, y+(dt/2)*k_2, u);
k_4 = func(t+dt  , y+    dt*k_3, u);

y = y + (dt/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
function s_next = rk4(f, t, s, dt)
    k1 = f(t, s) * dt;
    k2 = f(t + dt/2, s + k1/2) * dt;
    k3 = f(t + dt/2, s + k2/2) * dt;
    k4 = f(t + dt/2, s + k3/2) * dt;
    s_next = s + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
end
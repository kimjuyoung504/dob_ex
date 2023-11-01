clear; close all; clc;
addpath(genpath("./DynSysSim"))
addpath(genpath("./MissileSimEnv"))
addpath(genpath("./src"))

%% Run simulator
t0 = 0; % Simulation starts after burn-out. (Do not modify)
dt = 0.0025;
tf = 6.00;
ts = t0:dt:tf;

oc = 20; %omega_c = bandwidth
op_g = 4;
b0 = -300;
r = deg2rad(10);

limit = deg2rad(30);
Td = t0 + 0.05;
order = 1;
adrc = ADRC2("adrc", zeros(order + 1, 1));
adrc.setParams(b0, oc, order, op_g, r, limit, Td);

X = zeros(2, length(ts));
z = zeros(2, length(ts));
u = zeros(1, length(ts));

for i = 1:length(ts)
    u(:, i) = adrc.controller(ts(i), z(:,i));
    z(:, i+1) =  rk4(@adrc.dynEqns, ts(i), z(:, i), [u(:, i); X(1, i)], dt);
    X(:, i+1) = RK4_nl(@fx_nl, ts(i), X(:, i), u(i), dt);
end


figure(1)
subplot(2,1,1)
plot(ts, rad2deg(z(1, 1:end-1)), ts, rad2deg(X(1, 1:end-1)), ts, rad2deg(r) * ones(size(ts)))
legend('estimate', 'real', 'reference')
title('AOA')

subplot(2,1,2)
plot(ts, X(2, 1:end-1))
title('Q')


figure(2)
plot(ts, rad2deg(u), 'LineWidth', 1.5),title('input'), grid on, box on, legend('U')
rmpath(genpath("./DynSysSim"))
rmpath(genpath("./MissileSimEnv"))
rmpath(genpath("./src"))

function y = rk4(func, t, y, u, dt)
k_1 = func(t     , y           , u);
k_2 = func(t+dt/2, y+(dt/2)*k_1, u);
k_3 = func(t+dt/2, y+(dt/2)*k_2, u);
k_4 = func(t+dt  , y+    dt*k_3, u);

y = y + (dt/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);

end

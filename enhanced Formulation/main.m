clc; clear; clear all;


dt = 0.01; % Time step
tf = 20;
tspan = 0:dt:tf;

X0 = [0; 0; 0];

% x_hist = [];
dhat1_hist = [];
d_hist = [];
% u_hist = [];
X = X0;

lamb_1 = 1;
lamb_2 = 10;
k1 = 10;
k2 = 30;
for t = tspan
if t >= 6
    d = 5;
else
    d = 0;
end
d_hist = horzcat(d_hist, d);
    x1 = X(1);
    x2 = X(2);
    z = X(3);
    dhat1 = z + lamb_1 * x2;
    u = (-k1 * x1 - k2 * x2 - x1 * x2 - dhat1)/(1 + sin(x1)^2);

    % time history
    % x_hist = horzcat(x_hist, X(1:2));
    dhat1_hist = horzcat(dhat1_hist, dhat1);
    % u_hist = horzcat(u_hist, u);

    X_next1 = RK4(@system, t, X, u, dt, lamb_1);

    X = X_next1;
end
X = X0;
dhat2_hist = [];
% u_hist = [];
% x_hist = [];
for t = tspan

    x1 = X(1);
    x2 = X(2);
    z = X(3);
    dhat2 = z + lamb_2 * x2;
    u = (-k1 * x1 - k2 * x2 - x1 * x2 - dhat2)/(1 + sin(x1)^2);

    % time history
    % x_hist = horzcat(x_hist, X(1:2));
    dhat2_hist = horzcat(dhat2_hist, dhat2);
    % u_hist = horzcat(u_hist, u);

    X_next2 = RK4(@system, t, X, u, dt, lamb_2);
    X = X_next2;
end

figure
plot(tspan, d_hist, tspan, dhat1_hist(1, :), tspan, dhat2_hist(1, :),'LineWidth', 1.2), grid on
xlabel('Time (sec)')
ylabel('disturbance estimation')
legend('d','{\lambda} = 1','{\lambda} = 10')



function dX = system(t, X, u, lamb)
x1 = X(1);
x2 = X(2);
z = X(3);

% system dynamics
if t >= 6
    d = 5;
else
    d = 0;
end

dx1 = x2;
dx2 = x1 * x2 + (1 + sin(x1)^2) * u + d;

% disturbance observer
% lamb = param.lamb;
dz = -lamb * z - [0 lamb] * ([0; 1] * lamb * x2 ...
    + [x2; x1 * x2] + [0; (1 + sin(x1)^2) * u]);

dX = [dx1; dx2; dz];
end

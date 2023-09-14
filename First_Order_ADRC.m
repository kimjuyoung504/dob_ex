clear all; close all; clc;
%%First-Order ADRC with a First-Order Process

Tf = 3;
dt = 0.0025; %
tspan = 0:dt:Tf;

%Nominal parameter
K = 1;
T = 1;
b0 = K/T;
T_settle = 1;
Kp = 4;
s_eso = -40;
l_1 = -2 * s_eso;
l_2 = s_eso ^ 2;

A = [0 1; 0 0];
B = [b0; 0];
C = [1 0];
L = [l_1; l_2];
A_e = A - L * C;

param.A = A; param.B = B;
param.A_e = A_e; param.L = L;
param.T = 1; param.b0 = b0;
param.Kp = Kp;param.K = K;
x = zeros(2, length(tspan)); %
u = zeros(1, length(tspan));
u0 = zeros(1, length(tspan));
u(1) = 4;
x_e = zeros(2, length(tspan));
y = zeros(1, length(tspan));
r = 1;

for k = 1 : Tf/dt
    t = tspan(k);
    u0(k) = Kp * (r - x_e(1, k));
    u(k) = (u0(k) - x_e(2, k))/b0;
    param.y_h = x_e(1, k);
    param.d = 0.01 * exp( -0.01 * t ) * ( ( sin(15*t) )^2 * cos(15*t) );
    y( k + 1) = RK4(@fy, t, y(k), u(k), dt, param);
    param.y = y(k+1);
    
    x_e(:, k + 1) = RK4(@fx_e, t, x_e(:, k), u(k), dt, param);

end

figure(1)
plot(tspan, u),grid on, box on,xlabel('Time,sec'), ylabel('u');

figure(2)
plot(tspan, y),grid on, box on,xlabel('Time,sec'), ylabel('y');


function dy = fy(~, x, u, param)
% dy = param.Kp * (1 - param.y_h);
dy = -1/param.T * x + 1/param.T * param.d + (param.K/param.T + 1) * u;
end

function dx_e = fx_e(t, x_e, u, param)
dx_e = param.A_e * x_e + param.B * u + param.L * param.y;
end
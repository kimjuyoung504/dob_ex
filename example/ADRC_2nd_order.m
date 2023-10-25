clear all; close all; clc;
%%First-Order ADRC with a First-Order Process

Tf = 15;
dt = 0.0025; %
tspan = 0:dt:Tf;

%parameter
K = 1;
T = 1;
D = 1;
T_settle = 5;
s_cl = - 6/T_settle;
K_p = s_cl^2;
K_d = -2 * s_cl;
s_e = 10 * s_cl;
l_1 = -3 * s_e;
l_2 = 3 * s_e^2;
l_3 = - s_e^3;
b0 = K / T^2;

%Nominal parameter
A = [0 1 0; 0 0 1; 0 0 0];
B = [0; b0; 0];
C = [1 0 0];
L = [l_1; l_2; l_3];
A_e = A - L * C;
param.A_e = A_e; param.B = B;
param.L = L;param.T = T;
param.D = D; param.K = K;
param.b0 = b0;

x = zeros(3, length(tspan)); %
u = zeros(1, length(tspan));
u0 = zeros(1, length(tspan));
u(1) = 4;
x_e = zeros(3, length(tspan));
y = zeros(2, length(tspan));
r = 1;
a = ADRC();

for k = 1 : Tf/dt
    t = tspan(k);
    u0(k) = K_p * (r - x_e(1, k)) - K_d * x_e(2, k);
    u(k) = (u0(k) - x_e(3, k))/b0;
    param.y_h = x_e(1, k);
    param.d = 0.0001 * exp( -0.01 * t ) * ( ( sin(15*t) )^2 * cos(15*t) );
    % param.d = 0;
    y(:, k + 1) = RK4_param(@fy, t, y(:, k), u(k), dt, param);
    param.y = y(1, k+1);
    
    x_e(:, k + 1) = RK4_param(@fx_e, t, x_e(:, k), u(k), dt, param);

end

figure(1)
plot(tspan, u),grid on, box on,xlabel('Time,sec'), ylabel('u');

figure(2)
plot(tspan, y(1, :)),grid on, box on,xlabel('Time,sec'), ylabel('y');

function dy = fy(~, x, u, param)
% dy = param.Kp * (1 - param.y_h);
dy = [0 1;-1/param.T^2 -2 * param.D/param.T;] * x + [0; 1/param.T^2] * param.d + [0; param.b0] * u;
end

function dx_e = fx_e(t, x_e, u, param)
dx_e = param.A_e * x_e + param.B * u + param.L * param.y;
end
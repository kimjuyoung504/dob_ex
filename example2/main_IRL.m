clc; clear; close all;
addpath(genpath("./DynSysSim"))
addpath(genpath("./MissileSimEnv"))
addpath(genpath("./src"))



Q = diag([100 1]);
R = 100;
param.Q = Q;
param.R = R;

oc = 50; %omega_c = bandwidth
op_g = 3;
b0 = -1000;
r = deg2rad(10);
g1 = -57.7590;
g = [b0; g1];


dt = 0.0025;
tf = 6.00;
ts = 0:dt:tf;
T = round(0.01/dt); % round(0.04/dt); % integral reinforcement time: T*dt
N_intervals = 3; % control gain update period = N_intervals*T*dt

limit = deg2rad(30);
Td = 0.05;
order = 2;
adrc = ADRCIRL("adrc", zeros(order + 1, 1));
adrc.setParams(b0, oc, order, op_g, r, limit, Td, Q, R);

X = zeros(2, length(ts));
z = zeros(order + 1, length(ts));
u = zeros(1, length(ts));
u0 = zeros(1, length(ts));
J = zeros(1, length(ts));
K = zeros(1, length(ts));

X0 = deg2rad([5; -5]);
X(:, 1) = X0;
% IRL parameter
Phi = @(x) [x(1)^2 x(1)*x(2) x(2)^2].';
Phi_grad = @(x) [2*x(1) 0; x(2) x(1); 0 2*x(2)];

W0 = zeros(3, 1);
P0 = [W0(1) W0(2)/2; W0(2)/2 W0(3)];
K0 = R\g.' * P0;

W = W0;




Phidata = [];
Vdata = [];
TDerrordata = [];

phi = zeros(length(Phi(X0)),length(ts));



for i = 1:length(ts)
    
    phi(:, i) =  Phi(X(:, i));
    
    % IRL
    if (mod(i-1, T) == 0 && i ~= 1)   
        V_next  = J(i) - J(i-T) + W.'*phi(:,i);
        Phidata = cat(2, Phidata, phi(:, i-T));
        
        Vdata = cat(2, Vdata, V_next);
        
        W_cand = pinv(Phidata.') * Vdata.';
        
        if mod(i-1, T*N_intervals) == 0
            TDerrordata = cat(2, TDerrordata, max(abs(Vdata(:) - Phidata.'*W_cand)));
            if norm(Vdata) > 1e-5
                W = W_cand;
                Phidata = [];
                Vdata = [];

                % Policy Improvement
                P = [W(1) W(2)/2; W(2)/2 W(3)];
                K(:, i) =-0.5* R\ W.' * Phi_grad(X(:,i));
%                 K(:, i) = R\g.' * P;
            end
        end
    end
%     u0(i) = -K(:, i).' * X(:,i);
    u0(i) = - K(:, i);
% 
    % 4th RK
    u(:, i) = adrc.controller(ts(i), z(:,i), u0(i));
    z(:, i+1) = RK4_nl(@adrc.dynEqns, ts(i), z(:, i), [u(:, i); X(1, i)], dt);
    X(:, i+1) = RK4_nl(@fx_nl, ts(i), X(:, i), u(i), dt);
    J(i+1) = RK4(@fJ, [], J(i), [X(:,i); u(i)], dt, param);
    
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

figure(3)
plot(ts, rad2deg(u0), 'LineWidth', 1.5),title('input'), grid on, box on, legend('U0')

figure(4)
plot(ts, K(1,:), ts, K(2,:), 'LineWidth', 1.5),title('input'), grid on, box on, legend('K')

% figure(5)
% plot(ts, J, 'LineWidth', 1.5),title('Performance'), grid on, box on, legend('K')

rmpath(genpath("./DynSysSim"))
rmpath(genpath("./MissileSimEnv"))
rmpath(genpath("./src"))


% dot J
function Jdot = fJ(~, J, Z, param)
Jdot = Z.'*blkdiag(param.Q, param.R)*Z;
end



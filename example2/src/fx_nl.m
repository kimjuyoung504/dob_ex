function xdot = fx_nl(~, X, u)
% Aerodynamic coefficients
% bias_ratio1 = 0.8;
% bias_ratio2 = 0.8;
% bias_ratio3 = 0.8;

bias_ratio1 = 1;
bias_ratio2 = 1;
bias_ratio3 = 1;

aa = 0.3 * bias_ratio1;
am = 40.44 * bias_ratio3;
an = 19.373 * bias_ratio2;
bm = -64.015 * bias_ratio1;
bn = -31.023 * bias_ratio3;
cm = 2.922 * bias_ratio2;
cn = -9.717 * bias_ratio1;
dm = -11.803 * bias_ratio3;
dn = -1.948 * bias_ratio2;
em = -1.719 * bias_ratio1;

% Physical Configuration
S_ref = 0.0409;
d_ref = 0.2286;
mass = 204.02;
Iyy = 247.439;

% Fin dynamic model parameters
% Atmospheric coefficients (Ttable 2 in Ref.)
tro_alt = 10972.8; % troposphere boundary [m] (36000.0 ft)
str_alt = 20116.8; % stratosphere boundary [m] (66000.0 ft)
rho0sl = 1.1455; % Air density at sea level [kg/m^3] (2.377e-3 lb-s^2/ft^2)
rho0tr = 0.3415; % Air density at troposphere boundary [kg/m^3] (0.7086e-3 lb-s^2/ft^2)
Krhosl = 1.1103e-4; % [1/m] (3.36174e-5 1/ft)
KKrhotr = 1.5760e-4; % [1/m] (4.80377e-5 1/ft)
asl = 340.2787; % Sonic speed at the sea level [m/s] (1116.4 ft/sec)
atr = 295.0769; % Sonic speed at the troposphere boundary [m/s] (968.1 ft/sec)
Ka = 1.0358e-2; % [1/m] (0.00410833 1/ft)

% linear model parameter
h = 10000; %[m]
% rho = 0.012637;
% Vs = 340;%m/s
% mach = 2;
mach = 3.5;

rho = density(tro_alt, h,rho0sl, Krhosl, KKrhotr, rho0tr);
Vs = sonicSpeed(h, tro_alt, asl, Ka, atr);

alpha = X(1);
q = X(2);
fin = u(1);

Ca = aa;
Cz = an * alpha^3 + bn * alpha * abs(alpha) + cn * (2.0 - mach/3.0) * alpha + dn * fin;
Cm = am * alpha^3 + bm * alpha * abs(alpha) + cm * (-7.0 + 8.0*mach/3.0) * alpha + em * q + dm * fin;

f_1 = rho * Vs * mach * S_ref / (2.0 * mass) * (Cz * cos(alpha) - Ca * sin(alpha)) + q;
f_2 = rho * Vs^2 * mach^2 * S_ref * d_ref / (2.0 * Iyy) * Cm;

xdot =  [f_1; f_2];

    function out = density(tro_alt, h, rho0sl, Krhosl, KKrhotr, rho0tr)
        if h <= tro_alt
            out = rho0sl * exp(-Krhosl * h);
        else
            out = rho0tr * exp(-KKrhotr*(h - tro_alt));
        end
    end

    function out = sonicSpeed(h, tro_alt, asl, Ka, atr)
        if h <= tro_alt
            out =  asl - Ka * h;
        else
            out =  atr;
        end
    end

% function out = temperature(alt, h_trop, T0, L)
%     alt = min(max(alt, 0.0), h_trop);
%     out =  T0 - L * alt;
% end

end
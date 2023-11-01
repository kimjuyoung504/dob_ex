classdef MissileCoreDyn < DynSystems
    properties
    end
    methods
        function obj = MissileCoreDyn(name, x0, logOn)
            obj = obj@DynSystems(name, x0, logOn);
        end
        
        function dsdt = dynEqns(~, ~, s, force, moment, inertia, inertia_rate, mass)
            [u, v, w, p, q, r, ~] = disperse(s);
    
            a_x = force(1)/mass;
            a_y = force(2)/mass;
            a_z = force(3)/mass;

            M_x = moment(1);
            M_y = moment(2);
            M_z = moment(3);
            
            I_xx = inertia(1);
            I_yy = inertia(2);
            I_zz = inertia(3);
            I_xx_dot= inertia_rate(1);
            I_yy_dot= inertia_rate(2);
            I_zz_dot= inertia_rate(3);
            
            dudt = r * v - q * w + a_x;
            dvdt = p * w - r * u + a_y;
            dwdt = q * u - p * v + a_z;

            dpdt = M_x/I_xx - I_xx_dot/I_xx*p;
            dqdt = (I_zz - I_xx)/I_yy * p * r + M_y/I_yy - I_yy_dot/I_yy*q;
            drdt = (I_xx - I_yy)/I_zz * p * q + M_z/I_zz - I_zz_dot/I_zz*r;
            
            dmudt = p;
            dsdt = [dudt; dvdt; dwdt; dpdt; dqdt; drdt; dmudt];
        end
        
        function out = s2m(~, state)
            [u, v, w, p, q, r, mu] = disperse(state);
            V = sqrt(u^2+v^2+w^2);
            alpha = atan2(w, u);
            beta = asin(v/V);
            out = [alpha; beta; mu; p; q; r];
        end

        function [alpha, beta] = getIncideAngles(~, state)
            [u, v, w] = disperse(state(1:3));
            V = vecnorm(state(1:3), 2);
            alpha = atan2(w, u);
            beta = asin(v/V);
        end
    end
end
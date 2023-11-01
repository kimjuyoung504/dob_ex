% [Ref] S. Lee, Y. Kim, Y. Lee, G. Moon and B. -E. Jeon, "Robust-Backstepping Missile Autopilot Design Considering Time-Varying Parameters and Uncertainty," IEEE Transactions on Aerospace and Electronic Systems
% , Vol. 56, No. 6, pp. 4269-4287, Dec. 2020, doi: 10.1109/TAES.2020.2990819.
classdef ABSCMCtrl < DynSystems
    properties (Hidden)
        aerodyn
        k1
        k2 
        c1 
        Gamma1 
        Gamma2
        a 
        ks
        delta_lim
        alpha_lim = deg2rad(88.5) % In order to Prevent singularity issue of the inverse of g1
    end
    methods
        function obj = ABSCMCtrl(name, x0, logOn)
            obj = obj@DynSystems(name, x0, logOn);
        end
        
        function setParams(obj, aerodyn, k1, k2, c1, ks, a, Gamma1, Gamma2, delta_lim)
            obj.aerodyn = aerodyn;
            obj.k1 = k1;
            obj.k2 = k2;
            obj.c1 = c1;
            obj.ks = ks;
            obj.a = a;
            obj.Gamma1 = Gamma1;
            obj.Gamma2 = Gamma2;
            obj.delta_lim = delta_lim;
        end

        function dsdt = dynEqns(obj, e1, e2)
            Delta1_hat_dot = obj.Gamma1*e1(:);
            Delta2_hat_dot = obj.Gamma2*obj.c1.'*e2(:);
            dsdt = [Delta1_hat_dot; Delta2_hat_dot; e1(:); e2(:)];
        end
        
        function [u, x2d] = get_tvc_command(obj, t, x1, x2, ...
                Mach, angles_d, angles_d_dot, ...
                prop_data, air_data)

            [alpha, beta, ~] = disperse(x1);
            angles = x1;
            pqr = x2;
            Delta2 = obj.state(4:6);
            int_e1 = obj.state(7:9);
            int_e2 = obj.state(10:12);

            rho = air_data{1};
            Vs = air_data{2};
            S = air_data{3};
            d = air_data{4};
            V = Mach * Vs;
            Q = 1/2*rho*V^2;

            mass = prop_data{1};
            inertia = prop_data{2};
            inertia_rate = prop_data{3};
            Xcg = prop_data{4};
            thrust = prop_data{5};
            len_v = prop_data{6};


            % Aerodynamic coefficients
            delta = [0; 0; 0];
            fn_aerocoeff = @(alpha, beta, Mach, delta, t) obj.aerodyn.get_aerocoeff2(t, alpha, beta, Mach, delta, Xcg, d);
            aerocoeff = obj.aerodyn.get_aerocoeff2(t, alpha, beta, Mach, delta, Xcg, d);
            [Cx, Cy, Cz, Cl0, Clp, Cm0, Cmq, Cn0, Cnr, Cldp, ~, ~, ~, ~] = disperse(aerocoeff);
            x1 = angles;
            z1 = angles(:)-angles_d(:);
            ax = 1/mass * (thrust + Q*S*Cx);
            ay = 1/mass * Q*S*Cy;
            az = 1/mass * Q*S*Cz;
            f1 = obj.get_f1(alpha, beta, Mach, Vs, ax, ay, az);
            g1 = obj.get_g1(alpha, beta);
            
            dg1inv_dt = obj.get_dg1invdt(x1, x2, f1, g1);
            g1inv = obj.get_g1_inv(alpha, beta);

            z1dot = f1 + g1 * x2 - angles_d_dot;
            
            % Get x2d
            x2d = obj.get_x2_d(x1, Mach, angles_d,...
                    angles_d_dot, delta, Vs, Q, S, t,...
                    fn_aerocoeff, thrust, mass);
            
            % Inner loop design
            z2 = x2 - x2d;
            f2 = obj.get_f2(pqr, Q, inertia, inertia_rate, V, S, d, [Cl0, Cm0, Cn0, Clp, Cmq, Cnr]);
            len_t = Xcg - len_v;
            g2 = obj.get_g2_tvc(inertia, Q, S, d, len_t, thrust, Cldp);
            dx2ddt =  - dg1inv_dt * (f1 + obj.k1 * z1) ...
                      - g1inv * obj.k1 * z1dot;
            u = g2\( ...
                    -f2 ...
                    + dx2ddt ...
                    - obj.k2 * z2 ...
                    - obj.c1 * g1.'*z1...
                    - obj.ks * tanh(obj.a*z2) ...
                    - Delta2(:) ...
                ); 
            for i = 1:3 
                u(i) = max(min(u(i), obj.delta_lim(i)),-obj.delta_lim(i));
            end
        end
        
        function out = get_x2_d(obj, x1, x_int, angles_d, ...
                angles_d_dot, delta, Vs, Q, S, t, ...
                fn_aerocoeff, thrust, mass)
                
                z1 = x1(:)-angles_d(:);
                x_1_d_dot = angles_d_dot(:);
                alpha = x1(1);
                beta = x1(2);
                Mach = x_int(1);
                Delta1 = obj.state(1:3);

                aerocoeff = fn_aerocoeff(alpha, beta, Mach, delta, t);
                Cx = aerocoeff(1);
                Cy = aerocoeff(2);
                Cz = aerocoeff(3);
                ax = 1/mass*(thrust + Q*S*Cx);
                ay = 1/mass*(Q*S*Cy);
                az = 1/mass*(Q*S*Cz);
                f1 = obj.get_f1(alpha, beta, Mach, Vs, ax, ay, az);
                g1_inv = obj.get_g1_inv(alpha, beta);
                
                % Get x2d
                out = g1_inv*( ...
                            - f1 ...
                            + x_1_d_dot ...
                            - obj.k1*z1 ...
                            - Delta1(:) ...
                            );
        end

        function out = get_g1_inv(obj, alpha, beta)
            alpha = min(max(alpha, -obj.alpha_lim), obj.alpha_lim);
            out = [0, 0, 1;
                   1, -tan(alpha)*tan(beta), cos(alpha)*tan(beta) + sin(alpha)*tan(alpha)*tan(beta);
                   0, -1/cos(alpha), tan(alpha)];
        end
        
        function out = get_x2_d_dot(obj)

        end
    end
    methods (Static)
        function out = get_f1(alpha, beta, Mach, Vs, ax, ay, az)
            out = 1/(Mach*Vs)*[1/cos(beta)*(az*cos(alpha) - ax*sin(alpha));
                               -(ax*cos(alpha)*sin(beta)-ay*cos(beta) +az*sin(alpha)*sin(beta));
                               0];
        end

        function out = get_f2(pqr, Q, inertia, inertia_rate, V, S, d, aerocoeff)
            [p,q,r] = disperse(pqr);
            [I_xx, I_yy, I_zz] = disperse(inertia);
            [I_xx_dot, I_yy_dot, I_zz_dot] = disperse(inertia_rate);
            [Cl, Cm, Cn, Cl_p, Cm_q, Cn_r] = disperse(aerocoeff);
            out = [ Q*S*d/I_xx*(Cl + Cl_p*d/2/V*p) - I_xx_dot/I_xx*p;
                    (I_zz - I_xx)/I_yy*p*r + Q*S*d/I_yy*(Cm + Cm_q*d/2/V*q) - I_yy_dot/I_yy*q;
                    (I_xx - I_yy)/I_zz*p*q + Q*S*d/I_zz*(Cn + Cn_r*d/2/V*r) - I_zz_dot/I_zz*r];
        end

        function out = get_g1(alpha, beta)
            out = [-cos(alpha)*tan(beta), 1, -sin(alpha)*tan(beta);
                   sin(alpha),            0, -cos(alpha);
                   1,                     0, 0];
        end

        function out = get_g2(Cl_del, Cm_del_pitch, Cn_del_yaw, inertia, Q, S, d)
            [I_xx, I_yy, I_zz] = disperse(inertia);
            [Cl_del_roll, Cl_del_pitch, Cl_del_yaw] = disperse(Cl_del);
            out = Q*S*d*[Cl_del_roll/I_xx, Cl_del_pitch/I_xx, Cl_del_yaw/I_xx;
                         0,                Cm_del_pitch/I_yy, 0;
                         0,                0,                 Cn_del_yaw/I_zz];
        end

        function out = get_g2_tvc(inertia, Q, S, d, len_t, thrust, Cl_del_roll)
            thrust = max(thrust, 30000);
            [I_xx, I_yy, I_zz] = disperse(inertia);
            out = diag([Q*S*d*Cl_del_roll/I_xx; len_t*thrust/I_yy;-len_t*thrust/I_zz]);
        end

        function out = get_f_int(alpha, beta, Vs, ax, ay, az) % TODO: need to check
            out = 1/Vs*(ax * cos(alpha) * cos(beta) ...
                        + ay * sin(beta) ...
                        + az * sin(alpha) * cos(beta));
        end

        function out = get_dg1invdt(x1, x2, f1, g1)
            alpha = x1(1);
            beta = x1(2);
            dg11dx1 = zeros(3,3);
            dg12dx1 = [0,                      0,                            0;
                       -tan(beta)/cos(alpha)^2, -tan(alpha)/cos(beta)^2,     0;
                       -sin(alpha)/cos(alpha)^2,  0,                         0];

            dg13dx1 = [0,                       0, 0;
                       tan(alpha)/cos(alpha)*tan(beta), (cos(alpha) + sin(alpha)*tan(alpha))/cos(beta)^2, 0;
                       1/cos(alpha)^2, 0, 0];
            x1_dot = f1 + g1*x2;

            out = [dg11dx1*x1_dot, dg12dx1*x1_dot, dg13dx1*x1_dot];
        end
            
        function out = get_df1dx1(obj, alpha, beta, ax, ay, az)
           
        end

        function out = get_dfidx1(obj, alpha, beta, Mach, delta, Xcg, d)

            
        end

        function out = get_dDelta1dt(obj)
            
        end
        
        function out = da_dalpha(obj, alpha_t, beta, Mach, delta, Xcg, d, Q, S, mass)
            dCxyz_dab = obj.get_aerocoeffderiv(alpha_t, beta, Mach, delta, Xcg, d);
            dCxyz_da = dCxyz_dab(1, :).';
            out = Q*S/mass*dCxyz_da;
        end

        function out = da_dbeta(obj, alpha, beta_t, Mach, delta, Xcg, d, Q, S, mass)
            dCxyz_dab = obj.get_aerocoeffderiv(alpha, beta_t, Mach, delta, Xcg, d);
            dCxyz_db = dCxyz_dab(2, :).';
            out = Q*S/mass*dCxyz_db;
        end


        function out = get_dx2ddt(obj, x1, x1d, x1d_dot, x2, Mach, Vs, ax, ay, az)
            alpha = x1(1);
            beta = x1(2);
            Delta1 = obj.state(1:3);
            e1 = x1 - x1d;
            
            f1 = obj.get_f1(alpha, beta, Mach, Vs, ax, ay, az);
            g1 = obj.get_g1(alpha, beta);
            fi = obj.get_f_int(alpha, beta, mu, Vs, Mach, ax, ay, az);
            g1inv = obj.get_g1_inv(alpha, beta);
            dg1invdt = get_dg1invdt(x1, x2, Mach, Vs, ax, ay, az);
            df1dx1 = obj.get_df1dx();
            dfidxi = obj.get_dfidxi();
            Delta1_dot = obj.Gamma1*(f1 + g1*x2 - x1d_dot);

            out = dg1invdt*(-f1 + x1d_dot - obj.k1*e1 - Delta1(:)) +...
                  g1inv*(-df1dx1*(f1 + g1*x2) - dfidxi*fi - obj.k1*(f1 + g1*x2 - x1d_dot) - Delta1_dot);
        end

        function out = get_Phi1(alpha, beta, V, Q, S, Cx, mass, ax,ay,az)
            eta = [Q*S*Cx/mass, 0,  0,  -ax;
                   0,           ay, 0,  -ay;
                   0,           0,  az, -az];
            out = 1/V*[-sin(alpha)/cos(beta), 0, cos(alpha)/cos(beta);
                       -cos(alpha)*sin(beta), cos(beta), -sin(alpha)*sin(beta);
                       0, 0, 0]*eta;
        end
 
        function out = get_Phi2(pqr, aerocoeff, Q, S, d, inertia)
            [p,q,r] = disperse(pqr);
            [CL0, CM0, CN0] = disperse(aerocoeff);
            [I_xx, I_yy, I_zz] = disperse(inertia);
             out = [-Q*S*d/I_xx*CL0, 0, 0, Q*S*d/I_xx*CL0, 0, 0;
                    -I_xx/I_yy*p*r, - (I_zz-I_xx)/I_yy*p*r-Q*S*d/I_yy*CM0, I_zz/I_yy*p*r, 0, Q*S*d/I_yy*CM0, 0;
                     I_xx/I_zz*p*q, -I_yy/I_zz*p*q, -(I_xx-I_yy)/I_zz*p*q-Q*S*d/I_zz*CN0, 0, 0, Q*S*d/I_zz*CN0];
        end
        
        function out = get_Phi_int(alpha, beta, mu, Mach, Vs,Q,S,Cx,mass, ax, ay, az)
            eta = [Q*S*Cx/mass, 0,  0,  -ax;
                   0,           ay, 0,  -ay;
                   0,           0,  az, -az];
            out = 1/Vs*[cos(alpha)*cos(beta), sin(beta), sin(alpha)*cos(beta);
                        1/Mach*(sin(alpha)*cos(mu)+ cos(alpha)*sin(beta)*sin(mu)), -1/Mach*cos(beta)*sin(mu), -1/Mach*(cos(alpha)*cos(mu)-sin(alpha)*sin(beta)*sin(mu))]*eta;
        end
        
        function out = get_aerocoeffderiv(alpha0, beta0, Mach, delta, Xcg, d)
            out = obj.jac_fd(fn_aerocoeff, [alpha0; beta0]);
            function Cxyz = fn_aerocoeff(x)
                coeffs = obj.aerodyn.get_aerocoeff2(t, x(1), x(2), Mach, delta, Xcg, d);
                Cxyz = coeffs(1:3);
            end
        end

        function out = norm(x)
            out = vecnorm(x,2);
        end

        function jac = jac_fd(f, x, epsilon)
            % Calculate Jacobian of function f at given x
            % Standard finite difference method
            %
            % Inputs:
            %   f can be a vector of function, but make sure it is a row vector
            %   x is where the jacobian is being evaluated, it a row or column vector
            %   epsilon is a very small number

            if nargin < 3
                epsilon = 1e-5;
            end

            epsilon_inv = 1/epsilon;

            nx = length(x); % Dimension of the input x;

            f0 = feval(f, x); % caclulate f0, when no perturbation happens

            jac = zeros(length(f0), nx);

            % Do perturbation
            for i = 1 : nx
                xplus = x;
                xplus(i) =  x(i) + epsilon;
                jac(:, i) = (feval(f, xplus) - f0) .* epsilon_inv;
            end
        end
    end
end
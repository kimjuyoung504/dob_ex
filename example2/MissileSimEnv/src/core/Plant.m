classdef Plant < DynSystems
    properties (Hidden)
        dyn_pos
        dyn_uvwpqr
        dyn_quat
        aeEffects
        aerodyn
        airEnv
        prop
        lv = 2.71
        dataNames = {["\mu (deg)", "\gamma (deg)"],...
                     ["ax (m/s^2)", "ay (m/s^2)", "az (m/s^2)"], ...
                     ["\zeta_p (deg)", "\zeta_r (deg)"]}
        acc
    end
    
    methods
        function obj = Plant(name, s0, logOn, aerodatapath, propdatapath)   
            % Extract states
            NEDpos0 = s0(1:3);
            uvw0 = s0(4:6);
            pqr0 = s0(7:9);
            int_p0 = s0(10);
            quat0 = s0(11:14);

            % Set supporting subsystems
            airEnv = EarthAirEnv();
            aerodyn = Aerodynamic_SRAAM(aerodatapath);
            prop = Propulsion(propdatapath);
            
            % Set aeroelsticity
            modes = 1:5;

            % Set dynamical systems
            dyn_pos = Position("pos", NEDpos0, false);
            dyn_uvwpqr = MissileCoreDyn("coreDyn", [uvw0;pqr0;int_p0], false);
            dyn_quat = Quaternion("Quaternion", quat0, false);
            aeEffect_pitch = AeroelasticEffects("Aeroelasticity",modes, true);
            aeEffect_yaw = AeroelasticEffects("Aeroelasticity",modes, true);
            obj = obj@DynSystems(name, {dyn_pos, dyn_uvwpqr, dyn_quat, aeEffect_pitch, aeEffect_yaw}, logOn);

            obj.dyn_pos = dyn_pos;
            obj.dyn_uvwpqr = dyn_uvwpqr;
            obj.dyn_quat = dyn_quat;
            obj.aeEffects = {aeEffect_pitch, aeEffect_yaw};
            obj.aerodyn = aerodyn;
            obj.airEnv = airEnv;
            obj.prop = prop;
        end

        function obj = setParams(obj, prop_biases, aerocoeff_biases, noisestds, ignore)
            obj.prop.setParams(prop_biases);
            obj.aerodyn.setParams(aerocoeff_biases, noisestds);
            E = 0.4*1e9;
            d = 0.1294;
            Iyy = pi/64*d^4;
            Izz = pi/64*d^4;
            L = obj.lv;
            rho = obj.prop.getMass(0)/obj.lv; % kg /m^3
            area = 0.0132; % m^2
            xi = 0.0001;
            
            obj.aeEffects{1}.setParams(E, Izz, L, rho, area, xi, 0, ignore);
            obj.aeEffects{2}.setParams(E, Iyy, L, rho, area, xi, 0, ignore);
        end
        
        function dsdt = dynEqns(obj, t, s, input)
            % Get states
            [pos, uvwpqrmu, quat, ae_p, ae_y] = obj.splitStates(s);
            uvw = uvwpqrmu(1:3);
            Vm = vecnorm(uvwpqrmu(1:3),2);
            pqr = uvwpqrmu(4:6);
            alt = - pos(3);

            % Density
            rho = obj.airEnv.density(alt);
            Vs = obj.airEnv.soundSpeed(alt);

            % get proplusional parameters
            inertia = obj.prop.getInertia(t);
            inertia_rate = obj.prop.getInertiaRate(t);
            mass = obj.prop.getMass(t);
            Xcg = obj.prop.getXcg(t);
            
            % Velocity vector represented in NED coordinate
            dcmLocal2Body = obj.dyn_quat.quat2dcm(quat);
            dcmBody2Local = dcmLocal2Body.';
            NEDvel = dcmBody2Local * uvw;
            
            % Aerodynamic forces & moments
            delta_pqr = input;
            mach = Vm/Vs;
            mprop = t < obj.prop.toff;
            [alpha, beta] = obj.getIncidenceAngles(uvw);
            [q_pitch, q_pitch_dot] = obj.aeEffects{1}.get_q(t, ae_p);
            [q_yaw, q_yaw_dot] = obj.aeEffects{2}.get_q(t, ae_y);
            alpha_w = obj.aeEffects{1}.inducedIncidesAngles(uvw(1), alpha, Xcg, q_pitch, q_pitch_dot);
            beta_w = obj.aeEffects{2}.inducedIncidesAngles(uvw(1), beta, Xcg, q_yaw, q_yaw_dot);
            alpha_t = acos(cos(alpha_w) * cos(beta_w));
            phi_t = atan2(tan(beta_w), sin(alpha_w));
            [F_aero, M_aero] = obj.aerodyn.aeroForceAndMoment(alpha_t, phi_t, ...
                pqr, delta_pqr, Vm, Xcg, mprop, mach, rho);
            
            % Gravitational forces
            [phi, theta, ~] = obj.dyn_quat.quat2euler(quat);
            grav = obj.airEnv.grav();
            F_grav = grav*[-sin(theta)
                           sin(phi)*cos(theta)
                           cos(phi)*cos(theta)];
            
            % Thrust force & moment
            thrust = obj.prop.getThrustMag(t);
            delta_te = -obj.aeEffects{1}.dzeta_dl(obj.lv, q_pitch);
            delta_tr = obj.aeEffects{2}.dzeta_dl(obj.lv, q_yaw);
            F_thrust = thrust * [cos(delta_tr)*cos(delta_te);
                                 -sin(delta_tr);
                                 -cos(delta_tr)*sin(delta_te)];

            momentArm = obj.lv - Xcg;
            zetaT_p = obj.aeEffects{1}.elasticDisplacement(obj.lv, q_pitch);
            zetaT_y = obj.aeEffects{2}.elasticDisplacement(obj.lv, q_yaw);
            M_thrust = [0;
                        momentArm * thrust * sin(delta_te) * cos(delta_tr);
                        -momentArm * thrust * sin(delta_te)] ...
                      -[0;
                        F_thrust(3) * zetaT_p;
                        F_thrust(2) * zetaT_y];

            % Total force & moment
            force = F_aero + F_thrust + F_grav;
            moment = M_aero + M_thrust;
            
            % Generalized forces and moments for aeroelastic effects
            Fz = force(3);
            Fy = force(2);
            aeEffect_pitch = obj.aeEffects{1};
            aeEffect_yaw = obj.aeEffects{2};
            Gfy = aeEffect_yaw.generalizedForce(Fy);
            Gfz = aeEffect_pitch.generalizedForce(Fz);
            Msz = aeEffect_pitch.generalizedMass(mass);
            Msy = aeEffect_yaw.generalizedMass(mass);
            Gay = Gfy./Msy;
            Gaz = Gfz./Msz;
            zetaXcg_p = aeEffect_pitch.elasticDisplacement(Xcg, q_pitch);
            zetaXcg_y = aeEffect_yaw.elasticDisplacement(Xcg, q_yaw);
            
            % get body acceleration
            obj.acc = force/mass;

            % Vehicle Dynamical system Equations
            dpos_dt = obj.dyn_pos.dynEqns(t, NEDvel);
            duvwpqr_dt = obj.dyn_uvwpqr.dynEqns(t, uvwpqrmu, ...
                force, moment, inertia, inertia_rate, mass);
            dquat_dt = obj.dyn_quat.dynEqns(t, quat, pqr);
            dqs_pitch_dt = aeEffect_pitch.dynEqns(t, ae_p, Gaz);
            dqs_yaw_dt = aeEffect_yaw.dynEqns(t, ae_y, Gay);
            dsdt = [dpos_dt; duvwpqr_dt; dquat_dt; dqs_pitch_dt; dqs_yaw_dt];
            if obj.logData
                obj.data.fpa = rad2deg(obj.getFlightPathAngles(alpha, beta, phi, theta));
                obj.data.acc = obj.getAccBody(force, mass);
                obj.data.incides = rad2deg([alpha, beta]);
                obj.data.totalangles = rad2deg([alpha_t; phi_t]);
                obj.data.zetaXcg = rad2deg([zetaXcg_p, zetaXcg_y]);
                obj.data.zetadot = [aeEffect_pitch.data.dzetadt, aeEffect_yaw.data.dzetadt];
                obj.data.q_pitch = aeEffect_pitch.data.q;
                obj.data.dzetadl = [aeEffect_pitch.data.dzetadl, aeEffect_yaw.data.dzetadl];
                obj.data.Gay = Gay;
                obj.data.Gaz = Gaz;
                obj.data.incides_ae = rad2deg([alpha_w, beta_w]);
                obj.data.delta_thrust = rad2deg([delta_te; delta_tr]);
                obj.data.zetaT = rad2deg([zetaT_p, zetaT_y]);
            end
        end              

        function [alpha_t, phi_t] = getTotalAngles(obj, uvw)
            u = uvw(1);
            v = uvw(2);
            w = uvw(3);
            V = vecnorm(uvw, 2);
            alpha_t = acos(u/V);
            phi_t = atan2(v, w);
        end
        
        function [alpha, beta] = getIncidenceAngles(obj, uvw)
            V = vecnorm(uvw, 2);
            alpha = atan2(uvw(3), uvw(1));
            beta = asin(uvw(2)/V);
        end

        function [mu, gamma] = getFlightPathAngles(~, alpha, beta, phi, theta)
            sin_mu_cos_gamma = ( cos(theta)*sin(phi) )*cos(beta) + ( cos(alpha)*sin(theta)-sin(alpha)*cos(theta)*cos(phi) )*sin(beta);
            cos_mu_cos_gamma = ( cos(alpha)*cos(theta)*cos(phi) + sin(alpha)*sin(theta));
            mu = atan2(sin_mu_cos_gamma, cos_mu_cos_gamma);
            cos_gamma = cos_mu_cos_gamma/cos(mu);
            sin_gamma = ( cos(alpha)*sin(theta)-sin(alpha)*cos(theta)*cos(phi) )*cos(beta) - (cos(theta)*sin(phi))*sin(beta);
            gamma = atan2(sin_gamma, cos_gamma);
            if abs(mu) > pi/2
                mu = sign(mu)*pi - mu;
                gamma = sign(gamma)*pi - gamma;
            end
        end
    
        function out = getAccBody(~, force, mass)
            out = force/mass;
        end
    end

    methods (Static)
        function test()
            
            % Initial operating region
            alt0 = 11000;
            mach = 3.0;

            % Set initial state vector
            airEnv = EarthAirEnv();
            Vs = airEnv.soundSpeed(alt0);
            Vm = Vs * mach;
            alpha = deg2rad(0.0);
            beta = deg2rad(0.0);
            NEDpos0 = [0;0; -alt0];
            uvw0 = Vm * [ cos(alpha) * cos(beta);
                          sin(beta)
                          cos(alpha) * sin(beta)];
            quat0 = eul2quat([0,0,0], "ZYX");
            pqr0 = [0;0;0];
            int_p0 = 0;
            s0 = {NEDpos0, uvw0, pqr0, int_p0, quat0};

            % set dynamical system
            logOn = true;
            sys = Plant("dynsys", s0, logOn);
            
            % Control inputs
            delta_pqr = deg2rad([5;5;0]);

            % one step
            dt = 0.0025;
            s_next = sys.step(delta_pqr, dt);
        end
    end
end
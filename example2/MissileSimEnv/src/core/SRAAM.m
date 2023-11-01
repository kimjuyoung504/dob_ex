classdef SRAAM < DynSystems
    properties (Hidden)
        plant
        act
        dataNames = {["North (km)", "East (km)", "Down (km)"], ...
                     ["u (m/s)", "v (m/s)", "w (m/s)"], ...
                     ["p (deg/s)", "q (deg/s)", "r (deg/s)"], ...
                     ["\alpha_t (deg)", "\phi_t (deg)"], ...
                     ["ax (m/s^2)", "ay (m/s^2)", "az (m/s^2)"], ...
                     ["\zeta_p", "\zeta_y"],...
                     ["\zeta_p", "\zeta_y"],...
                     ["q_1","q_2","q_3","q_4","q_5"],...
                     ["\zeta_p", "\zeta_y"],...
                     ["Gay1","Gay2","Gay3","Gay4","Gay5"],...
                     ["Gaz1","Gaz2","Gaz3","Gaz4","Gaz5"]}
    end

    methods
        function obj = SRAAM(name, s0, logOn, aerodatapath, propdatapath)
            if nargin == 2
                logOn = false;
            end
            plant0 = s0(1:14);
            plant = Plant("plant", plant0, true, aerodatapath, propdatapath);
            act0 = s0(15:22);
            tailfin = FourTailFin("actuator", act0, false);
            subsystems = {plant, tailfin};
            obj = obj@DynSystems(name, subsystems, logOn);
            obj.plant = plant;
            obj.act = tailfin;
        end

        function obj = setParams(obj, prop_biases, aerocoeff_biases, noisestds, ignore)
            if nargin == 1
                prop_biases = zeros(9,1);
                aerocoeff_biases = 0;
                noisestds = 0;
                ignore = true;
            end
            obj.plant.setParams(prop_biases, aerocoeff_biases, noisestds, ignore);
        end

        function dsdt = dynEqns(obj, t, s, input)
            % Actuator commands
            delta_pqr_cmd = input;  

            % Obtain derivatives of subsystems
            [s_plant, s_act] = obj.splitStates(s);
            delta_pqr = obj.act.getDelta_pqr(s_act([1,3,5,7]));
            d_plant = obj.plant.dynEqns(t, s_plant, delta_pqr);
            d_act = obj.act.dynEqns(t, s_act, delta_pqr_cmd);
            dsdt = [d_plant; d_act];
            if obj.logData
                m2km = 0.001;
                r2d = 180/pi;
                obj.data.pos = s(1:3)*m2km;
                obj.data.velBody = s(4:6);
                obj.data.pqr = s(7:9) * r2d;
                obj.data.totalangles = obj.plant.data.totalangles;
                obj.data.acc = obj.plant.data.acc;
                obj.data.zeta = obj.plant.data.zetaXcg;
                obj.data.zetadot = obj.plant.data.zetadot;
                obj.data.q = obj.plant.data.q_pitch;
                obj.data.dzetadl = obj.plant.data.dzetadl;
                obj.data.Gay = obj.plant.data.Gay;
                obj.data.Gaz = obj.plant.data.Gaz;
            end
        end
    
        function out = getAcc(obj, t, s)
            % Get states
            [s_plant, s_act] = obj.splitStates(s);
            [pos, uvwpqrmu, quat, ae_p, ae_y] = obj.plant.splitStates(s_plant);
            uvw = uvwpqrmu(1:3);
            Vm = vecnorm(uvwpqrmu(1:3),2);
            pqr = uvwpqrmu(4:6);
            alt = - pos(3);

            % Density
            rho = obj.plant.airEnv.density(alt);
            Vs = obj.plant.airEnv.soundSpeed(alt);

            % get proplusional parameters
            mass = obj.plant.prop.getMass(t);
            Xcg = obj.plant.prop.getXcg(t);
            
            % Aerodynamic forces & moments
            delta_pqr = obj.act.fin2ctrl(s_act([1,3,5,7]));
            mach = Vm/Vs;
            mprop = t < obj.plant.prop.toff;
            [alpha, beta] = obj.plant.getIncidenceAngles(uvw);
            [q_pitch, q_pitch_dot] = obj.plant.aeEffects{1}.get_q(t, ae_p);
            [q_yaw, q_yaw_dot] = obj.plant.aeEffects{2}.get_q(t, ae_y);
            alpha_w = obj.plant.aeEffects{1}.inducedIncidesAngles(uvw(1), alpha, Xcg, q_pitch, q_pitch_dot);
            beta_w = obj.plant.aeEffects{2}.inducedIncidesAngles(uvw(1), beta, Xcg, q_yaw, q_yaw_dot);
            alpha_t = acos(cos(alpha_w) * cos(beta_w));
            phi_t = atan2(tan(beta_w), sin(alpha_w));
            [F_aero, ~] = obj.plant.aerodyn.aeroForceAndMoment(alpha_t, phi_t, ...
                pqr, delta_pqr, Vm, Xcg, mprop, mach, rho);
            
            % Gravitational forces
            [phi, theta, ~] = obj.plant.dyn_quat.quat2euler(quat);
            grav = obj.plant.airEnv.grav();
            F_grav = grav*[-sin(theta)
                           sin(phi)*cos(theta)
                           cos(phi)*cos(theta)];
            
            % Thrust force & moment
            thrust = obj.plant.prop.getThrustMag(t);
            delta_te = -obj.plant.aeEffects{1}.dzeta_dl(obj.plant.lv, q_pitch);
            delta_tr = obj.plant.aeEffects{2}.dzeta_dl(obj.plant.lv, q_yaw);
            F_thrust = thrust * [cos(delta_tr)*cos(delta_te);
                                 -sin(delta_tr);
                                 -cos(delta_tr)*sin(delta_te)];

            % Total force & moment
            force = F_aero + F_thrust + F_grav;
            out = force/mass;
        end
    end

    methods (Static)
        function test()
            % Initial operating region
            alt0 = 15000;
            mach = 3.0;

            % Set initial state vector
            airEnv = EarthAirEnv();
            Vs = airEnv.soundSpeed(alt0);
            Vm = Vs * mach;
            alpha = deg2rad(0.0);
            beta = deg2rad(0.0);
            NEDpos0 = [0;0; -alt0];
            uvw0 = Vm * [ cos(alpha) * cos(beta);
                          sin(beta);
                          sin(alpha) * cos(beta)];
            quat0 = eul2quat([0,0,0], "ZYX").';
            pqr0 = [0;0;0];
            int_p0 = 0;
            s0 = [NEDpos0; uvw0; pqr0; int_p0; quat0; zeros(8,1)];

            % Set dynamical system
            logOn = true;
            aerodataPath = "./Data/aeroDB_SRAAM.mat";
            propdataPath = "./Data/propulsionDB.mat";
            sys = SRAAM("dynsys", s0, logOn, aerodataPath, propdataPath);
            sys.setParams();

            % Control inputs
            delta_pqr_cmd = deg2rad([0;10;0]);

            % Propagation
            t0 = 3.0;
            dt = 0.0025;            
            tf = 6.0;
            tspan = t0:dt:tf;
            sim = Simulator(sys).propagate(tspan, delta_pqr_cmd);
            sim.report();
            log = sim.log;
            
            log.acc.subplots(1, "Acceleration")
            log.pos.subplots(1,"Position");
            log.velBody.subplots(2, "Velocity");
            log.pqr.subplots(3, "Angular rate");
            log.totalangles.subplots(4, "Total angles");
            log.zeta.subplots(5);
            log.zetadot.subplots(6);
            log.q.subplots(7);
            log.dzetadl.subplots(8);
            log.Gay.plot(9);
            log.Gaz.plot(10);
        end
    end
end
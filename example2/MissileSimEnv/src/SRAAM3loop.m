classdef SRAAM3loop < DynSystems
    properties (Hidden)
        sraam
        ctrl
        dataNames = {["a_y (m/s^2)", "a_z (m/s^2)"],...
                     ["a_y (m/s^2)", "a_z (m/s^2)"],...
                     ["\alpha (deg)", "\beta (deg)"],...
                     ["\alpha (deg)", "\beta (deg)"],...
                     ["\zeta_p (m)" , "\zeta_y (m)"],...
                     ["\zeta_p /dl", "\zeta_y /dl"],...
                     ["\delta_{Te} (deg)","\delta_{Tr} (deg)"],...
                     ["\zeta_p (m)", "\zeta_y (m)"]}
    end

    methods
        function obj = SRAAM3loop(name, s0, logOn)
            aerodatapath = "./core/Data/aeroDB_SRAAM.mat";
            propdatapath = "./core/Data/propulsionDB.mat";
            sraam = SRAAM("sraam", s0, logOn, aerodatapath, propdatapath);
            ctrl = ThreeLoopAutopilot("rollPitchYaw", [0;0], logOn);
            obj = obj@DynSystems(name, {sraam, ctrl}, logOn);
            obj.sraam = sraam;
            obj.ctrl = ctrl;
        end
        
        function obj = setParams(obj, prop_bias, aerocoeff_bias, noiseStd, ignore)
            obj.sraam.setParams(prop_bias, aerocoeff_bias, noiseStd, ignore);
            deltaCmdLim = deg2rad(30);
            filenameThrustOff = "./core/Data/ThreeLoopGain_FuelOff.mat";
            filenameThrustOn = "./core/Data/ThreeLoopGain_FuelOn.mat";
            obj.ctrl.setParams(deltaCmdLim, filenameThrustOff, filenameThrustOn);

        end

        function dsdt = dynEqns(obj, t, s, u)
            accCmds = u;

            [s_sraam, s_ctrl] = obj.splitStates(s);
            alt = -s_sraam(3);
            p = s_sraam(7);
            q = s_sraam(8);
            r = s_sraam(9);
            pqr = [p;q;r];
            uvw = s_sraam(4:6);
            quat = s_sraam(11:14);
            e_p = s_ctrl(1);
            e_y = s_ctrl(2);
            eulZYX = quat2eul(quat(:).');
            phi = eulZYX(3);
            Vs = obj.sraam.plant.airEnv.soundSpeed(alt);
            mach = vecnorm(uvw, 2)/Vs;
            
            % 3-loop controller
            accMax = obj.ctrl.getGainFromTable(t, mach, alt);
            delta_p_cmd = obj.ctrl.rollController(phi, p);
            delta_q_cmd = obj.ctrl.pitchController(e_p, q);
            delta_r_cmd = obj.ctrl.yawController(e_y, r);
            delta_pqr_cmd = [delta_p_cmd; delta_q_cmd; delta_r_cmd];
            d_sraam = obj.sraam.dynEqns(t, s_sraam, delta_pqr_cmd);
            ay = obj.sraam.plant.acc(2); % [m/s^2]
            az = obj.sraam.plant.acc(3); % [m/s^2]
            d_ctrl = obj.ctrl.dynEqns(t, s_ctrl, accCmds, ay, az, pqr, accMax);
            dsdt = [d_sraam; d_ctrl];
            if obj.logData
                obj.data.acc = obj.sraam.data.acc(2:3);
                obj.data.acc_cmd = flip(u*9.81);
                obj.data.incides = obj.sraam.plant.data.incides;
                obj.data.incides_ae = obj.sraam.plant.data.incides_ae;
                obj.data.zeta_Xcg = obj.sraam.plant.data.zetaXcg;
                obj.data.dzetadl = obj.sraam.data.dzetadl;
                obj.data.delta_thrust = obj.sraam.plant.data.delta_thrust;
                obj.data.zeta_T = obj.sraam.plant.data.zetaT;
            end
        end
    end

    methods (Static)
        function test()
%             set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.5,'defaultfigurewindowstyle','docked');
%             set(0,'defaultfigurecolor','white');
%             set(0, 'DefaultLineLineWidth', 1);
            addpath(genpath("./"))
            % Initial operating region
            alt0 = 4000;
            mach = 1.0;

            % Set initial state vector
            airEnv = EarthAirEnv();
            Vs = airEnv.soundSpeed(alt0);
            Vm = Vs * mach;
            alpha = deg2rad(0.0);
            beta = deg2rad(0.0);
            NEDpos0 = [0;0; -alt0];
            uvw0 = Vm * [ cos(alpha) * cos(beta);
                          sin(beta);
                          cos(alpha) * sin(beta)];
            quat0 = eul2quat([0,0,0], "ZYX").';
            pqr0 = [0;0;0];
            int_p0 = 0;
            s0 = [NEDpos0; uvw0; pqr0; int_p0; quat0; zeros(8,1)];

            % Uncertainties
            prop_bias = [0;0;0;0;0;0;0;0;0];
            aerocoeff_bias = 0;
            noiseStd = 0;
            
            % Aeroelastic effects
            aeOff = false;

            % Set systems
            logOn = true;
            sys = SRAAM3loop("sys", s0, logOn);
            sys.setParams(prop_bias, aerocoeff_bias, noiseStd, aeOff);

            % Run simulation
            ts = 0:0.0025:5;
            input = @(t) [sys.stepCmd(t, 1, 5),...
                          sys.stepCmd(t, 1, 5)]; % [grav]
            sim = Simulator(sys).propagate(ts, input);
            sim.report();
            log = sim.log;
            

            log.acc.subplots(1, "Acceleration histories"); 
            plt = log.acc_cmd.subplots(1);
            for i = 1:length(plt)
                plt{i}.Color= "k";
                plt{i}.LineStyle = "--";
            end
%             log.incides.subplots(2);
            log.incides_ae.subplots(2, "Incidence angles");
            log.zeta_Xcg.subplots(3, "Elastic displacement at the center of gravity");
            log.zeta_T.subplots(4, "Elastic displacement at the nozzle")
            log.dzetadl.subplots(5, "d\zeta /dl");
            log.delta_thrust.subplots(6,"TVC nozzle deflection");
            
            rmpath(genpath("./"))
        end

        function cmd = stepCmd(t, stepTimes, stepValues)
            cmd = 0;
            if nargin <= 1
                return
            end

            idx = sum(t >= stepTimes);
            if idx == 0
                return
            end
            cmd = stepValues(idx);
        end
    end
end
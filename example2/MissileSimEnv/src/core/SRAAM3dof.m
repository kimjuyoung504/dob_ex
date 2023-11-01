classdef SRAAM3dof < SRAAM
    methods
        function obj = SRAAM3dof(name, s0, logOn, aerdatapath, propdatapath)
            [u, w, q, theta, h] = disperse(s0);
            quat = eul2quat([0 theta 0]).';
            NEDpos = [0;0;-h];
            uvw = [u;0;w];
            pqr = [0; q; 0];
            fins = zeros(8,1);
            s0 = [NEDpos; uvw; pqr; 0; quat; fins];
            obj = obj@SRAAM(name, s0, logOn, aerdatapath, propdatapath);
        end

        function dsdt = dynEqns(obj, t, s, del_q_cmd)
            del_pqr_cmd = [0; del_q_cmd; 0];
            dsdt = obj.dynEqns@SRAAM(t, s, del_pqr_cmd);
            if obj.logData
                obj.data.acc = obj.plant.data.acc;
            end
        end

        function out = fullstates(obj, s_long)
            [u, w, q, theta, alt, dr, del_q, fin_dot] = disperse(s_long);
            NEDpos = [dr;0;-alt];
            uvw = [u;0;w];
            pqr = [0;q;0];
            fins = obj.act.ctrl2fin([0;del_q;0]);
            acts = nan(8,1);
            for i = 1:4
                acts(2*i-1) = fins(i);
                acts(2*i) = fin_dot;
            end
            quat = eul2quat([0 theta 0]).';
            ae = zeros(2*obj.plant.aeEffects{1}.n_modes, 1);
            state = [NEDpos; uvw; pqr; 0; quat; ae; ae; acts];
            out = state;
        end

        function out = reducedstates(obj, s)
            dr = s(1);
            alt = -s(3);
            u = s(4);
            w = s(6);
            q = s(8);
            quat = s(11:14);
            eulZYX = quat2eul(quat.');
            theta = eulZYX(2);
            fins = s([35,37,39,41]);
            delta_pqr = obj.act.fin2ctrl(fins);
            del_q = delta_pqr(2);
            fin_dot = s(36);
            out = [u, w, q, theta, alt, dr, del_q, fin_dot].';
        end
        
        function az = getAcc(obj, t, s)
            s = obj.fullstates(s);
            acc = obj.getAcc@SRAAM(t, s);
            az = acc(3);
        end

        function s_next = step_rk4(obj, t, s, input, dt)
            s = obj.fullstates(s);
            s_next = obj.step(t, s, input, dt);
            s_next = obj.reducedstates(s_next);
        end
    end

    methods (Static)
        function test()
            addpath("./Data/")
            aerdatapath = "./Data/aeroDB_SRAAM.mat";
            propdatapath = "./Data/propulsionDB.mat";
            logOn = true;

            % initial states
            alt = 15000;
            mach = 3.0;
            airEnv = EarthAirEnv();
            Vs = airEnv.soundSpeed(alt);
            Vm = Vs * mach;
            alpha = deg2rad(0.0);
            u = Vm*cos(alpha);
            w = Vm*sin(alpha);
            q = 0;
            theta = 0;
            downrange = 0;
            del_q = 0; % Fin deflection angle of longitudinal plane
            fin_dot = 0;
            s0 = [u; w; q; theta; alt; downrange; del_q; fin_dot];

            % Set dynamical system
            sys = SRAAM3dof("dyn", s0, logOn, aerdatapath, propdatapath);
            sys.setParams();
            sim = Simulator(sys);

            t0 = 3.0;
            dt = 0.0025;
            tf = 6.0;
            ts = t0:dt:tf;
            input = deg2rad(10);
            sim.propagate(ts, input);
            sim.report();
            log = sim.log;

            log.pos.subplots(1,"Position");
            log.velBody.subplots(2, "Velocity");
            log.pqr.subplots(3, "Angular rate");
            log.totalangles.subplots(4, "Total angles");
            log.acc.subplots(5, "Acceleration");
            rmpath("./Data/")
        end
    end
end
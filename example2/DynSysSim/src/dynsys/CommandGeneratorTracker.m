classdef CommandGeneratorTracker < DynSystems
    properties (Hidden)
        A
        B
        dataNames = ["x", "xdot"]
    end

    methods
        function obj = CommandGeneratorTracker(name, x0, logOn)
            obj = obj@DynSystems(name, x0, logOn);            
        end

        function obj = setParams(obj, A, B)
            obj.A = A;
            obj.B = B;
        end

        function dsdt = dynEqns(obj, ~, s, u)
            dsdt = obj.A * s + obj.B * u;
            if obj.logData
                obj.data.state = s;
            end
        end
    end

    methods (Static)
        function test()
            addpath("../simulator/")
            addpath("../utils/")
            addpath("../solvers/")
            addpath("../datalogger/")

            mag = 2;
            x0 = [mag; 0];
            omega = 2;
            A = [0 1;-omega^2 0];
            B = [0; 0];
            logOn = true;
            sys = CommandGeneratorTracker("test", x0, logOn).setParams(A, B);
            dt = 0.01;
            tf = 10.0;
            sim = Simulator(sys).propagate(dt, tf, @(t) 0.0);
            sim.report();
            log = sim.log;

            log.state.subplots(1, "state"); hold on;
            subplot(211);
            tspan = 0:dt:tf;
            sol = @(t) mag*cos(omega*t);
            plot(tspan, sol(tspan),"k--"); 
            legend("IVP solution", "Analytic solution");
            subplot(212);
            sol = @(t) -mag*omega*sin(omega*t);
            plot(tspan, sol(tspan), "k--");
            legend("IVP solution", "Analytic solution");

            rmpath("../simulator/")
            rmpath("../utils/")
            rmpath("../solvers/")
            rmpath("../datalogger/")
        end
    end
end
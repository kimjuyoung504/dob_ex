% This is second-order system class which is originated from source code of Prof. Shin
classdef NonlinearSecondOrderSystem < SecondOrderSystem
%     properties (Hidden)
%         dataNames = {["x","xdot"],["xdot","xddot"],"Command"}
%     end    
        
    methods
        function obj = NonlinearSecondOrderSystem(name, s0, logOn)  
            obj = obj@SecondOrderSystem(name, s0, logOn);
            obj.name = name;
        end
        
        function sDot = dsdt(obj, ~, s, u) 
            t = obj.time;
            x = s(1);
            x_dot = s(2);
            x_ddot = - obj.sat(x, obj.xLim) * obj.omega^2 ...
                     - 2 * obj.xi*obj.omega*obj.sat(x_dot, obj.xDotLim)...
                     + u * obj.omega^2;
            sDot = [x_dot; x_ddot];
            if obj.checkLoggerOn(t)
                obj.data.state = s;
                obj.data.stateDot = sDot;
                obj.data.cmd = u;
            end
        end                     
    end
    
    methods (Static)
        function test()
            % add path
            addpath("../simulator/")
            addpath("../utils/")
            addpath("../solvers/")
            addpath("../datalogger/")

            xi = 0.7;
            omega = 50; % [rad/s]
            xLim =  deg2rad(30); %[rad]
            xDotLim = deg2rad(500); % [rad/s]
            s0 = [0;0];
            logOn = true;
            sys = NonlinearSecondOrderSystem('sys', s0, logOn).setParams(xi, omega, xLim, xDotLim);   
            dt = 0.0005;
            tf = 2;
            input = @(t) stepCmd(t, 0, deg2rad(10));

            sim = Simulator(sys).propagate(dt, tf, input);
            sim.report();
            log = sim.log;

            % analytic solution
            c = @(t, xi, omega, v) v*(1 - exp(-xi*omega*t)/sqrt(1-xi^2)...
                .*sin( omega*sqrt(1-xi^2)*t + atan( sqrt(1-xi^2)/xi)));
            traj = c(sim.tspan, xi, omega, deg2rad(10));
             
            % Plots
            log.cmd.plot(1); hold on;
            log.state.plot(1, 'Second order system history', 1);
            plot(sim.tspan, traj,"--")
            legend("Command","Ours","Analytic solutions")
            
            % Remove path
            rmpath("../simulator/")
            rmpath("../utils/")
            rmpath("../solvers/")
            rmpath("../datalogger/")
        end
        function out = sat(x, lim)
            out = min(max(x, -lim),lim);
        end
    end
end
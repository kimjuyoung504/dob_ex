classdef SecondOrderSystem < DynSystems
    properties (Hidden)
        xi
        omega
        xLim = inf 
        xDotLim = inf      
        dataNames = {["x","xdot"],["xdot","xddot"],"Command"}
    end    
        
    methods
        function obj = SecondOrderSystem(name, s0, logOn)            
                obj = obj@DynSystems(name, s0, logOn);
        end     
        
        function obj = setParams(obj, xi, omega, xLim, xDotLim)
            obj.xi = xi;
            obj.omega = omega; % [rad/s]
            obj.xLim = xLim;
            obj.xDotLim = xDotLim;
        end

        function sDot = dynEqns(obj, ~, s, u) 
            s1 = s(1);
            s2 = s(2);
            A = [0, 1; -obj.omega^2, -2*obj.xi*obj.omega];
            B = [0; obj.omega^2];
            sDot = A * s(:) + B* u;
            
            if abs(s2) >= obj.xDotLim
                si = sign(s2);
                sDot(2) = si * min(si * sDot(2), 0);
                obj.state(2) = si * obj.xDotLim;
                sDot(1) = obj.state(2);
            end
            
            if abs(s1) >= obj.xLim
                si = sign(s1);
                obj.state(1) = si * obj.xLim;
                sDot(1) = si * min(si * sDot(1), 0);
                obj.state(2) = sDot(1);
                if sDot(1) == 0
                    sDot(2) = si * min(si * sDot(2), 0);
                end               
            end
            
            if obj.logData
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
            omega = 250;
            xLim =  deg2rad(10); %[rad]
            xDotLim = deg2rad(500); % [rad/s]
            logOn = true;
            sys = SecondOrderSystem('sys', [0;0], logOn).setParams(xi, omega, xLim, xDotLim);
            dt = 0.0025;
            tf = 2;
            input = @(t) stepCmd(t, [0,0.5,1], deg2rad([5,-5,5]));
            
            sim = Simulator(sys).propagate(dt, tf, input);
            sim.report();
            log = sim.log;

            log.state.plot(1, 'Second order system history', 1);
            log.cmd.plot();

            % remove path
            rmpath("../simulator/")
            rmpath("../utils/")
            rmpath("../solvers/")
            rmpath("../datalogger/")
        end
    end
end
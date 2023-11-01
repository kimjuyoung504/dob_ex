% [Ref.1] K. A. Wise, and D. J. B roy, “Agile missile Dynamics and Control,” Journal of Guidance, Control, and Dynamics, Vol.21, No. 3, pp.441-449, 1998.
%
classdef ThrustVector < DynSystems
    properties
        dataNames = {{'delta elevator (deg)','delta rudder (deg)'},...
                {'delta rate elevator (deg/s)', 'detla rate rudder (deg/s)'}}
    end
    methods
        function obj = ThrustVector(name, s0, logOn)
            names = {'TVC_fin_elevator', 'TVC_fin_rudder'};            
            nozzles = cell(1,2);
            for i = 1:2
                nozzles{i} = SecondOrderSystem(names{i}, s0(2*i-1:2*i), logOn);
            end
            obj = obj@DynSystems(name, nozzles, logOn);
        end
        
        function obj = setParams(obj, xi, omega, delta_lim, delta_rate_lim)
            for i = 1:obj.subSysNum
               obj.subSysCell{i}.setParams(xi, omega, delta_lim, delta_rate_lim);
            end
        end

        function sDot = dynEqns(obj, t, s, u)
            if nargin == 1
                u = obj.input(obj.time);
            end
            sDot = nan(4, 1);
            for i = 1:2
                sDot(2*i-1:2*i, 1) = obj.subSysCell{i}.dynEqns(t, s(2*i-1:2*i), u(i));
            end
            if obj.logData
                obj.data.delta = rad2deg(obj.state([1,3]));
                obj.data.delta_rate = rad2deg(obj.state([2,4]));
            end
        end

        function [force, moment] = actuate(obj, thrust, momentArm)
            assert(isscalar(thrust), 'Thrust must be scalar')
            assert(thrust >=0, 'Thrust must be non-negative')

            % refer to [Ref.1]
            de = obj.state(1);
            dr = obj.state(3);
            force = thrust * [cos(de)*cos(dr);
                              -sin(dr);
                              -sin(de)*cos(dr)];
            moment = [0;
                      momentArm*thrust*sin(de)*cos(dr);
                      -momentArm*thrust*sin(dr)];
        end

        function [de,dr] = get_angles(obj)
            de = obj.state(1);
            dr = obj.state(3);
        end
    end

    methods (Static)
        function test()
            addpath('../utils')
            
            logOn = true;
            sys = ThrustVector('tvc', zeros(4,1), logOn).setParams(0.7, 100, deg2rad(10), deg2rad(100));
            dt = 0.0025;
            tf = 1.0;
            v = deg2rad([5,-5,5]);
            
            input = @(t) [stepCmd(t, [0, 0.2,0.4], v), stepCmd(t, [0,0.2,0.4], v)];
            
            sim = Simulator(sys).propagate(dt, tf, input);
            sim.report();
            log = sim.log;

            % plot
            log.delta.subplots(1,'Nozzle deflection angles');
            log.delta_rate.subplots(2, "Nozzle deflection angle rates");
            
            rmpath('../utils');
        end
    end
end

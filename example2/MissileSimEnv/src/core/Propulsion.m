 classdef Propulsion < handle
    properties (Hidden)   
        preEvaluate = false
        thrustMap
        IxxMap
        IyyMap
        IzzMap
        IxxdotMap
        IyydotMap
        IzzdotMap
        massMap
        XcgMap
        thTbl
        IxxTbl
        IyyTbl
        IzzTbl
        IxxdotTbl
        IyydotTbl
        IzzdotTbl
        massTbl
        XcgTbl
        thrustBias = 0
        IxxBias = 0
        IyyBias = 0
        IzzBias = 0
        IxxdotBias = 0
        IyydotBias = 0
        IzzdotBias = 0
        massBias = 0
        XcgBias = 0
    end

    properties (Constant)
        toff = 2.69
    end
    
    methods
        function obj = Propulsion(propulsionDB, timespan)            
            load(propulsionDB,'time','thrust','I_xx','I_zz', ...
                'I_xx_dot', 'I_zz_dot', 'Mass', 'Xcg');
            time = time(:);
            obj.thTbl = DataTable(thrust(:), time);  
            obj.IxxTbl = DataTable(I_xx(:), time);
            obj.IyyTbl = DataTable(I_zz(:), time);
            obj.IzzTbl = DataTable(I_zz(:), time);
            obj.IxxdotTbl = DataTable(I_xx_dot(:), time);
            obj.IyydotTbl = DataTable(I_zz_dot(:), time);
            obj.IzzdotTbl = DataTable(I_zz_dot(:), time);
            obj.massTbl = DataTable(Mass(:), time);
            obj.XcgTbl = DataTable(Xcg(:), time);

            if nargin >= 2
                ts = round(min(timespan, max(time)), 6);  
                obj.thrustMap = containers.Map(ts, obj.thTbl.interpolate(ts));
                obj.IxxMap = containers.Map(ts, obj.IxxTbl.interpolate(ts));
                obj.IyyMap = containers.Map(ts, obj.IyyTbl.interpolate(ts));
                obj.IzzMap = containers.Map(ts, obj.IzzTbl.interpolate(ts));
                obj.IxxdotMap = containers.Map(ts, obj.IxxdotTbl.interpolate(ts));
                obj.IyydotMap = containers.Map(ts, obj.IyydotTbl.interpolate(ts));
                obj.IzzdotMap = containers.Map(ts, obj.IzzdotTbl.interpolate(ts));
                obj.massMap = containers.Map(ts, obj.massTbl.interpolate(ts));
                obj.XcgMap = containers.Map(ts, obj.XcgTbl.interpolate(ts));
                obj.preEvaluate = true;
            end           
        end

        function obj = setParams(obj, biases)
            obj.thrustBias = biases(1);
            obj.IxxBias = biases(2);
            obj.IyyBias = biases(3);
            obj.IzzBias = biases(4);
            obj.IxxdotBias = biases(5);
            obj.IyydotBias = biases(6);
            obj.IzzdotBias = biases(7);
            obj.massBias = biases(8);
            obj.XcgBias = biases(9);
        end

        function out = getThrustMag(obj, t)
            if obj.preEvaluate == true
                t = round(t, 6);
                out =  obj.thrustMap(t);
            else
                out = obj.thTbl.interpolate(t);
            end
            out = out * ( 1 + obj.thrustBias );
        end

        function out = getInertia(obj, t)
            if obj.preEvaluate == true
                t = round(t, 6);
                Ixx = obj.IxxMap(t);
                Iyy = obj.IyyMap(t);
                Izz = obj.IzzMap(t);
            else
                Ixx = obj.IxxTbl.interpolate(t);
                Iyy = obj.IyyTbl.interpolate(t);
                Izz = obj.IzzTbl.interpolate(t);
            end           
            Ixx = Ixx * ( 1 + obj.IxxBias );
            Iyy = Iyy * ( 1 + obj.IyyBias );
            Izz = Izz * ( 1 + obj.IzzBias );
            out = [Ixx; Iyy; Izz];
        end

        function out = getInertiaRate(obj, t)
            if obj.preEvaluate == true
                t = round(t, 6);
                Ixxdot = obj.IxxdotMap(t);
                Iyydot = obj.IyydotMap(t);
                Izzdot = obj.IzzdotMap(t);
            else
                Ixxdot = obj.IxxdotTbl.interpolate(t);
                Iyydot = obj.IyydotTbl.interpolate(t);
                Izzdot = obj.IzzdotTbl.interpolate(t);
            end
            Ixxdot = Ixxdot * ( 1 + obj.IxxdotBias );
            Iyydot = Iyydot * ( 1 + obj.IyydotBias );
            Izzdot = Izzdot * ( 1 + obj.IzzdotBias );
            out = [Ixxdot, Iyydot, Izzdot];
        end

        function out = getMass(obj, t)
            if obj.preEvaluate == true
                t = round(t, 6);
                out =  obj.massMap(t);
            else
                out = obj.massTbl.interpolate(t);
            end 
            out = out * ( 1 + obj.massBias );
        end

        function out = getXcg(obj, t)
            if obj.preEvaluate == true
                t = round(t, 6);
                out =  obj.XcgMap(t);
            else
                out = obj.XcgTbl.interpolate(t);
            end    
            out = out * ( 1 + obj.XcgBias);
        end

        function out = getThrustMagData(obj, t)
            if obj.preEvaluate == true
                t = round(t, 6);
                out =  obj.thrustMap(t);
            else
                out = obj.thTbl.interpolate(t);
            end           
        end

        function [Ixx, Iyy, Izz] = getInertiaData(obj, t)
            if obj.preEvaluate == true
                t = round(t, 6);
                Ixx = obj.IxxMap(t);
                Iyy = obj.IyyMap(t);
                Izz = obj.IzzMap(t);
            else
                Ixx = obj.IxxTbl.interpolate(t);
                Iyy = obj.IyyTbl.interpolate(t);
                Izz = obj.IzzTbl.interpolate(t);
            end           
        end

        function [Ixxdot, Iyydot, Izzdot]= getInertiaRateData(obj, t)
            if obj.preEvaluate == true
                t = round(t, 6);
                Ixxdot = obj.IxxdotMap(t);
                Iyydot = obj.IyydotMap(t);
                Izzdot = obj.IzzdotMap(t);
            else
                Ixxdot = obj.IxxdotTbl.interpolate(t);
                Iyydot = obj.IyydotTbl.interpolate(t);
                Izzdot = obj.IzzdotTbl.interpolate(t);
            end
        end

        function out = getMassData(obj, t)
            if obj.preEvaluate == true
                t = round(t, 6);
                out =  obj.massMap(t);
            else
                out = obj.massTbl.interpolate(t);
            end           
        end

        function out = getXcgData(obj, t)
            if obj.preEvaluate == true
                t = round(t, 6);
                out =  obj.XcgMap(t);
            else
                out = obj.XcgTbl.interpolate(t);
            end           
        end
    end

    methods (Static)
        function test()
            PATH_propulsionDB = './Data/propulsionDB.mat';
            tspan = 0:0.01:6;
            p = Propulsion(PATH_propulsionDB, tspan);
            thrust = nan(length(tspan), 1);
            Ixx_dot = nan(length(tspan), 1);
            Iyy_dot = nan(length(tspan), 1);
            Izz_dot = nan(length(tspan), 1);
            for i = 1:length(tspan)
                thrust(i) = p.getThrustMag(tspan(i));
                [Ixx_dot(i), Iyy_dot(i), Izz_dot(i)] = p.getInertiaRate(tspan(i));
            end
            plot(tspan, thrust)
            plot(tspan, Ixx_dot, tspan, Iyy_dot, tspan, Izz_dot)
        end
    end
end
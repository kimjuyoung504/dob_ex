classdef Simulator < handle
    properties
        system
        log = []
        isStopCondSatisfied = false
    end
    properties (Hidden)
        INIT = 1
        elapsedTime = 0
        logSwitches
        logHist
        fieldNames
        tspan
    end
  
    methods
        function obj = Simulator(system, logHist)
            if nargin <= 1
                logHist = false;
            end
            classes = superclasses(system);
            if ~strcmp(classes{end-1}, 'DynSystems')
                error('invalid system')
            end
            obj.system = system;
            obj.logHist = logHist;
        end

        function obj = propagate(obj, tspan, u)
            if isempty(obj.system.state)
                error('state should be initialized before calling propagate method');
            end
            if isnumeric(u)
               u = @(t) u; 
            end

            obj.tspan = round(tspan, 6);
            numiter = length(tspan);
            dt = tspan(2) - tspan(1);
            obj.system.updateTimes(tspan(1));
            tic
            for i = 1:numiter
                t = obj.system.time();
                s = obj.system.state();
                s_next = obj.system.step(t, s, u(t), dt);
                t_next = round(t + dt ,6);  
                obj.stackdata();
                obj.system.updateState(s_next);
                obj.system.updateTimes(t_next);
            end
            obj.elapsedTime = toc;                    
            if obj.INIT 
                obj.INIT = false;
            end
        end

        function set_ZOHCountLimit(obj, n)
            obj.system.ZOHCountLimit = n;
            obj.system.ZOHCallCount = n-1;
        end

        function stackdata(obj)
            data = obj.system.data;
            if isempty(obj.log) && ~isempty(data)
                fnames = fieldnames(data);
                fieldLength = length(fnames);
                [~, obj.log] = DataInventory(fnames, obj.system.dataNames, obj.tspan);
                for i = 1:fieldLength
                    obj.log.(fnames{i}).append(data.(fnames{i}));
                end
                obj.fieldNames = fnames;
            elseif ~isempty(data)
                fnames = obj.fieldNames;
                for i = 1:length(fnames)
                    obj.log.(fnames{i}).append(data.(fnames{i}));
                end
            end
        end
        
        
        function report(obj)
            fprintf('[Elapsed simulation time] %.2f seconds \n',obj.elapsedTime)
        end
        
        function f = set_f(obj, u)
            if nargin == 1
                f = @(t, s) obj.system.stateSpaceEqn(t, s);
            elseif isa(u, 'function_handle')
                f = @(t, s) obj.system.stateSpaceEqn(t, s, u(t));
            else
                f = @(t, s) obj.system.stateSpaceEqn(t, s, u);
            end
        end
    end                
end
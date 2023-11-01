classdef DynSystems < handle
    properties   
        name
        time = 0
        state 
        data
    end

    properties (Hidden)
        stateDeriv
        stateSize
        subSysCell
        subSysNum = 0
        subSysSize 
        subSysIdxes
        isLogOn = false
        logData = false
        holder
    end
        
    methods
        function obj = DynSystems(name, in, logOn) 
            if nargin == 3
                obj.isLogOn = logOn;
            end
            assert(or(iscell(in), isnumeric(in)), 'invalid input');
            if iscell(in)              
                obj.subSysCell = in;
                obj.subSysNum = numel(in);  
                startIdx = 1;
                for i = 1:obj.subSysNum
                    stSize = in{i}.stateSize;
                    obj.subSysSize(i) = stSize;
                    stateArr(startIdx:startIdx + stSize -1, 1) = in{i}.state;
                    startIdx = startIdx  + stSize;
                end                   
                obj.updateState(stateArr);
                obj.stateSize = startIdx - 1;
            else
                obj.updateState(in);   
                obj.stateSize = numel(in);
            end
            obj.name = name;
        end

        function out = dynEqns(obj, t, ~, u)
            if obj.subSysNum == 0
                error("Create the 'dynEqns' method in the "+obj.name+" dynamical system")
            end
            out = nan(sum(obj.subSysSize), 1);
            startIdx = 1;
            for i = 1:obj.subSysNum
                subsys = obj.subSysCell{i};
                s = subsys.state;
                stSize = subsys.stateSize;
                if obj.logData == true
                    subsys.logData = true;
                    out(startIdx:startIdx + stSize-1) = subsys.dynEqns(t, s, u{i});
                    obj.data.(subsys.name) = subsys.data;
                    subsys.logData = false;
                else
                    out(startIdx:startIdx + stSize-1) = subsys.dynEqns(t, s, u{i});
                end
                startIdx = startIdx + stSize;
            end
        end

        function s_next = step(obj, t, s, u, dt)
            obj.holder = [];
            %t = obj.time;
            %s = obj.state;
            f = obj.setODEfun(u);
            
            obj.switchLogData(true);
            k1 = f(t, s) * dt;
            obj.switchLogData(false);
            k2 = f(t + dt/2, s + k1/2) * dt;
            k3 = f(t + dt/2, s + k2/2) * dt;
            k4 = f(t + dt, s + k3) * dt;
            s_next = s + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
        end

        function [varargout] = ZeroOrderHold(obj, varargin)
            if isempty(obj.holder)
                varargout = varargin;
                obj.holder = varargin;
            else
                varargout = obj.holder;
            end
        end

        function f = setODEfun(obj, u)
            if nargin == 1
                f = @(t, s) obj.dynEqns(t, s);
            else
                f =  @(t, s) obj.dynEqns(t, s, u);
            end
        end

        function [varargout] = getSubStates(obj)
            varargout = cell(obj.subSysNum, 1);
            for i = 1:obj.subSysNum
                varargout{i} = obj.subSysCell{i}.state;
            end
        end

        function updateState(obj, state)
            startIdx = 1;
            for i = 1: obj.subSysNum
                sysSize = obj.subSysSize(i);
                obj.subSysCell{i}.updateState(state(startIdx:startIdx + sysSize - 1));
                startIdx = startIdx + sysSize;
            end            
            obj.state = state;
        end
        
        function updateTimes(obj, t)
            obj.time = t;
            for i = 1:obj.subSysNum
                obj.subSysCell{i}.updateTimes(t);
            end
        end

        function reset(obj, state)
            obj.state = state;
            startIdx = 1;
            for i = 1:obj.subSysNum
                subSysLen = obj.subSysSize(i);
                subState = state(startIdx: startIdx + subSysLen - 1);
                obj.subSysCell{i}.reset(subState);
                startIdx = startIdx + subSysLen;
            end
            obj.time = 0;
            obj.data = [];
        end

        function [varargout] = splitStates(obj, state)
            stSizes = obj.subSysSize;
            varargout = cell(obj.subSysNum, 1);
            startIdx = 1;
            for i = 1:obj.subSysNum
                varargout{i} = state(startIdx:startIdx + stSizes(i) -1);
                startIdx = startIdx  + stSizes(i);
            end
        end

        function switchLogData(obj, bool)
            if obj.isLogOn == true
                obj.logData = bool;
                for i = 1:obj.subSysNum
                    subsys = obj.subSysCell{i};
                    if subsys.isLogOn == true
                        subsys.switchLogData(bool);
                    end
                end
            end
        end
    end
end
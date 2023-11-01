classdef MonteCarloSimulator < handle
    properties
        system
        sweepVars
        sweepVarNum
        simNum
        dataSet
        flattenVar
        log
    end
    methods
        function obj = MonteCarloSimulator(systemFcnHandle, varargin)
            obj.system = systemFcnHandle;
            obj.sweepVars = varargin(:);
            obj.sweepVarNum = numel(varargin);
            % flatten sweep variables
            obj.flatSweepVars;
        end
        
        function preallocate(obj, fieldNames, dt, tf)
            timespan = 0:dt:tf;
            dataLen = length(timespan);
            obj.log = MatDataInventory(fieldNames, dataLen, obj.simNum, timespan);
        end
        
        function run(obj, dt, tf, input, fieldNames)
            dataLen = length(0:dt:tf);
            obj.preallocate(fieldNames, dt, tf);
            
            progBar = ProgressBar(obj.simNum, 'UpdateRate', 100,'Title', 'Simulating');
            for i = 1:obj.simNum
                flatVar = num2cell(obj.flattenVar(i, :));
                sys = obj.system(flatVar{:});
                sys.setMonteCarloSim(dataLen, fieldNames)
                sim = Simulator(sys, 'rk4', [1 0]).propagate(dt, tf, input);
                obj.log.matData(:, i, :) = sim.system.MClog.data;
                progBar(1, [], []);
            end
            progBar.release();
            obj.log.mat2str;
        end
        
        function parRun(obj, dt, tf, input, sysNames, fieldNames)
            dataLen = round(tf/dt);
            obj.preallocate(sysNames, fieldNames, dataLen);

            % initialize systems
            sys = nan(obj.simNum, 1);
            for i = 1:obj.simNum
                flatVar = num2cell(flatVars(i, :));
                sys(i) = obj.system(flatVar{:});           
            end
            
%             parfor i = 1:obj.simNum
%                 sim = Simulator(sys(i), 'rk4', true).propagate(dt, tf, input);
%                 log = sim.History;
%                 for j = 1:dataNum
%                     dataset(:, i, j) = log.(sysNames{j}).(fieldNames{j}).data
%                 end
%             end
            obj.dataSet = dataset;
        end
        
        function flatSweepVars = flatSweepVars(obj)
            grid = cell(obj.sweepVarNum, 1);
            [grid{:}] = ndgrid(obj.sweepVars{:});
            obj.simNum = numel(grid{1});
            flatSweepVars = nan(obj.simNum, obj.sweepVarNum);
            for i = 1:obj.sweepVarNum
                flatSweepVars(:, i) = reshape(grid{i},[],1);
            end
            obj.flattenVar = flatSweepVars;
        end
    end
end
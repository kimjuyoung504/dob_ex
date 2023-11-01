classdef Position < DynSystems
    methods
        function obj = Position(name, s0, logOn)
            obj = obj@DynSystems(name, s0, logOn);
        end

        function dsdt = dynEqns(~, ~, vel)
            dsdt = vel;
        end
    end
end
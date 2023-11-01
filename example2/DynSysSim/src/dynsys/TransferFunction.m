classdef TransferFunction < DynSystems
    properties (Hidden)
        order
        denominatorCoeff
        numeratorCoeff
        dataNames
    end
    methods    
        function obj = TransferFunction(name, x0, logOn)
            obj = obj@DynSystems(name, x0, logOn);
        end
        
        function obj = setParams(obj, order, denomCoeff, numerCoeff)
            denomNum = numel(denomCoeff);
            numerNum = numel(numerCoeff);
            assert(and(denomNum == order + 1, numerNum <= order), 'invalid number of coeiffcients');
            if numerNum < order
                dif = order - numerNum;
                zeroCoeffs = zeros(1, dif);
                numerCoeff = [zeroCoeffs numerCoeff(:).'];
            end            
            obj.denominatorCoeff = denomCoeff(:).';
            obj.numeratorCoeff = numerCoeff(:).';
            obj.order = order;
        end

        function dsdt = dynEqns(obj, t, s, u)
            assert(length(obj.numeratorCoeff) == length(u),'the dimension of input dose not match with the dim. of state')
            u = u(:);
            an = obj.denominatorCoeff(1);
            yn = ( obj.numeratorCoeff * u(end:-1:1) - obj.denominatorCoeff(2:end) * s(end:-1:1) )/ an; 
            dsdt = [s(2:end); yn];
            if obj.logData
                obj.data.state = s;
                obj.data.stateDot = dsdt;
            end
        end

        function out = get.dataNames(obj)
            stateNames = cell(1, obj.order);
            stateDerivNames = stateNames;
            for i = 1:obj.order
                stateNames{i} = ['state',num2str(i)];
                stateDerivNames{i} = ['stateDot',num2str(i)];
            end
            out = {stateNames, stateDerivNames};
        end
    end
    
    methods (Static)
        function test()
            addpath("../simulator/")
            addpath("../utils/")
            addpath("../solvers/")
            addpath("../datalogger/")

            % case 1 : first order lag system
            tau = 0.1;
            order = 1;
            denominator = [tau 1];
            numerator = 1;
            logOn = true;
            sys = TransferFunction('firstOrderSys', 0, logOn).setParams(order, denominator, numerator);
            dt = 0.0025;
            tf = 5;
            input = @(t) stepCmd(t, 0, 1);
            sim1 = Simulator(sys).propagate(dt, tf, input);
            log = sim1.log;

            log.state.subplots(1, 'First order system state history');
            log.stateDot.subplots(2, 'First order system state derivative history');
            

            % case 2: second order system
            xi = 0.7;
            omega = 50; %[rad/s]
            order = 2;
            denominator = [1 2*xi*omega omega^2];
            numerator = omega^2;
            sys = TransferFunction('secondOrderSys', [0;0], logOn).setParams(order, denominator, numerator);
            dt = 0.0025;
            tf = 5;
            input = @(t) [stepCmd(t, 0, 1); 0];
            sim2 = Simulator(sys).propagate(dt, tf, input);
            log = sim2.log;
            
            log.state.subplots(3, 'second order system state history');
            log.stateDot.subplots(4, 'second order system state derivative history');
            
            % case 3: second order system with zero                        
            sys = TransferFunction('secondOrderSysWithZero', [0;0], logOn).setParams(order, [0.009 0.33 1], [-0.0363 1]);            
            sim3 = Simulator(sys).propagate(dt, tf, input);
            log = sim3.log;
            
            log.state.subplots(5, 'second order system with zero state history');
            log.stateDot.subplots(6, 'second order system with zero state derivative history');
            
            rmpath("../simulator/")
            rmpath("../utils/")
            rmpath("../solvers/")
            rmpath("../datalogger/")
        end
    end
end
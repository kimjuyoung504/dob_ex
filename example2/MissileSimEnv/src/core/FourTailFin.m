 classdef FourTailFin < DynSystems
    properties        
        dataNames = {["fin1", "fin2", "fin3", "fin4"],...
                     ["fin rate 1", "fin rate 2", "fin rate 3","fin rate 4"],...
                     ["\delta p", "\delta q", "\delta r"],...
                     ["\delta p", "\delta q", "\delta r"]}
    end
    
    methods
        function obj = FourTailFin(name, x0, logOn)
            names = {'fin1','fin2','fin3','fin4'};
            xi = 0.7;
            omega = 150;
            deltaLim = deg2rad(30);
            deltaRateLim = deg2rad(500);
            fins = cell(4, 1);
            for i = 1:4
                fins{i} = SecondOrderSystem(names{i}, x0(2*i-1:2*i), false);
                fins{i}.setParams(xi, omega, deltaLim, deltaRateLim);
            end            
            obj = obj@DynSystems(name, fins, logOn);
        end
        
        function sDot = dynEqns(obj, t, s, u)
            sDot = nan(8, 1);
            finCmd = obj.ctrl2fin(u);
            for i =1:4
                sDot(2*i-1:2*i, 1) = obj.subSysCell{i}.dynEqns(t, s(2*i-1:2*i), finCmd(i));
            end
            
            if obj.logData
                fin = s([1, 3, 5, 7]);
                obj.data.fin = rad2deg(fin);
                obj.data.finRate = rad2deg(s([2,4,6,8]));
                obj.data.delta_pqr = rad2deg(obj.getDelta_pqr(fin));
                obj.data.delta_pqr_cmd = rad2deg(u);
            end          
        end
         
        function out = getDelta_pqr(obj, fin)
            out = fin2ctrl(obj, fin);
        end
        
        function finCmd = ctrl2fin(obj, delta_pqr)
            finCmd = [ -1, 1, -1 ; 
                       -1, 1, 1 ;  
                       1, 1, -1 ; 
                       1, 1, 1 ] * delta_pqr(:);                   
        end
        
        function ctrlCmd = fin2ctrl(obj, finCmd)
            ctrlCmd = [ -0.2500   -0.2500    0.2500    0.2500 ;
                        0.2500    0.2500    0.2500    0.2500 ; 
                        -0.2500    0.2500   -0.2500    0.2500 ] * finCmd(:);
        end
    end
    
    methods (Static)
        function test()
            addpath("../")
            % set four tail fin system
            x0 = zeros(8,1);
            logOn = true;
            aeroFinSys = FourTailFin('fourTailFin', x0, logOn);            
            t0 = 0;
            dt = 0.0025;
            tf = 5;
            ts = t0:dt:tf;
            
            % set command
            input = @(t) stepCmd(t, [1,2], [1,-1]) * deg2rad([0, 1, 0]);
            sim = Simulator(aeroFinSys).propagate(ts, input);
            log = sim.log;

            % plot
            log.fin.subplots(1,'Fin deflection');
            log.delta_pqr.subplots(2, '\delta pqr'); hold on;
            log.delta_pqr_cmd.subplots(2);
            rmpath("../")
        end               
    end
end
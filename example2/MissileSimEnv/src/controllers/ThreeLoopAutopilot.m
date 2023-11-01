classdef ThreeLoopAutopilot < DynSystems
    properties (Hidden)        
        omega = 10
        zeta = 0.7
        pole = 20
        Ka = 0.06
        Kr =  0.001
        Wi = 0.01
        Kdc =  0.01
        rKx
        rKp
        errorPitch = 0
        errorYaw = 0
        shapedCmdDot = 0
        shapedCmd = 0
        deltaCmdLim = deg2rad(30)
        gainsThrustOn
        gainsThrustOff      
        gravConst = 9.81;
        altLim
        machLim
        dataNames = ['Kdc', 'Ka', 'Wi', 'Kr'];
    end      
    
    methods
        function obj = ThreeLoopAutopilot(name, x0, logOn)
            obj = obj@DynSystems(name, x0, logOn);
        end
        
        function obj = setParams(obj, deltaCmdLim, filenameThrustOff, filenameThrustOn)
            obj.deltaCmdLim = deltaCmdLim;
            if nargin >= 3
                load(filenameThrustOff, 'rKp', 'rKx', 'Ko', 'Ka', 'Wi', 'Kr', 'ALT', 'MACH', 'AccelMax');
                samplePoints = {MACH, ALT};
                obj.altLim = [min(ALT); max(ALT)];
                obj.machLim = [min(MACH); max(MACH)];
                obj.gainsThrustOff.Kdc= DataTable(Ko, samplePoints{:});
                obj.gainsThrustOff.Ka = DataTable(Ka, samplePoints{:});
                obj.gainsThrustOff.Wi = DataTable(Wi, samplePoints{:});
                obj.gainsThrustOff.Kr = DataTable(Kr, samplePoints{:});
                obj.gainsThrustOff.AccelMax = DataTable(AccelMax, samplePoints{:});
                obj.gainsThrustOff.rKp = DataTable(rKp, MACH, ALT);
                obj.gainsThrustOff.rKx = DataTable(rKx, MACH, ALT);

                if nargin >= 4
                    load(filenameThrustOn, 'rKp', 'rKx', 'Ko', 'Ka', 'Wi', 'Kr', 'AccelMax');                
                    obj.gainsThrustOn.Kdc= DataTable(Ko, samplePoints{:});
                    obj.gainsThrustOn.Ka = DataTable(Ka, samplePoints{:});
                    obj.gainsThrustOn.Wi = DataTable(Wi, samplePoints{:});
                    obj.gainsThrustOn.Kr = DataTable(Kr, samplePoints{:});
                    obj.gainsThrustOn.AccelMax = DataTable(AccelMax, samplePoints{:});
                    obj.gainsThrustOn.rKp = DataTable(rKp, MACH, ALT);
                    obj.gainsThrustOn.rKx = DataTable(rKx, MACH, ALT);
                end
            end
        end

        function dsdt = dynEqns(obj, ~, ~, u, ay, az, pqr, AccMax)
            acczCmd = u(1); % [g]
            accyCmd = u(2); % [g]
            [acczCmd, accyCmd] = obj.commandScaling(AccMax, [acczCmd, accyCmd]);

            q = pqr(2); % [rad/s]
            r = pqr(3); % [rad/s]

            errorPitchDot = obj.pitchDyn(az, q, acczCmd);
            errorYawDot = obj.yawDyn(ay, r, accyCmd);
            dsdt = [errorPitchDot; errorYawDot];

            if obj.logData
                obj.data.threeLoopGains = [obj.Kdc, obj.Ka, obj.Wi, obj.Kr];
            end
        end
        
        function [Kdc, Ka, Wi, Kr] = setGain(obj, alpha)
            % airframe linearized model
            transFcnCoff = obj.vehicle.linearizedDynamicModel(alpha);
            [Kdc, Ka, Wi, Kr]= obj.polePlacing(transFcnCoff);

            % outputs
            obj.Kdc = Kdc;
            obj.Ka = Ka;
            obj.Wi = Wi;
            obj.Kr = Kr;
        end

        function [AccMax, Kdc, Ka, Wi, Kr, rKx, rKp] = getGainFromTable(obj, t, mach, alt)
            if obj.altLim(1) > alt || obj.altLim(2) < alt
                warning('Current altitude must be within the [%s, %s]', num2str(obj.altLim(1)), num2str(obj.altLim(2)))
            end
            if obj.machLim(1) > mach || obj.machLim(2) < mach
                warning('Current mach number must be within the [%s, %s]', num2str(obj.machLim(1)), num2str(obj.machLim(2)))
            end
            mach = min(max(mach, obj.machLim(1)), obj.machLim(2));
            alt = min(max(alt, obj.altLim(1)), obj.altLim(2));
            if ~isempty(obj.gainsThrustOn)
                if t <= 2.69
                    Kdc = obj.gainsThrustOn.Kdc.interpolate(mach, alt);
                    Ka = obj.gainsThrustOn.Ka.interpolate(mach, alt);
                    Wi = obj.gainsThrustOn.Wi.interpolate(mach, alt);
                    Kr = obj.gainsThrustOn.Kr.interpolate(mach, alt);
                    AccMax = obj.gainsThrustOn.AccelMax.interpolate(mach, alt);
                    rKx = obj.gainsThrustOn.rKx.interpolate(mach, alt);
                    rKp = obj.gainsThrustOn.rKp.interpolate(mach, alt);
                else
                    Kdc = obj.gainsThrustOff.Kdc.interpolate(mach, alt);
                    Ka = obj.gainsThrustOff.Ka.interpolate(mach, alt);
                    Wi = obj.gainsThrustOff.Wi.interpolate(mach, alt);
                    Kr = obj.gainsThrustOff.Kr.interpolate(mach, alt);
                    AccMax = obj.gainsThrustOff.AccelMax.interpolate(mach, alt);
                    rKx = obj.gainsThrustOff.rKx.interpolate(mach, alt);
                    rKp = obj.gainsThrustOff.rKp.interpolate(mach, alt);
                end
            else
                Kdc = obj.gainsThrustOff.Kdc.interpolate(mach, alt);
                Ka = obj.gainsThrustOff.Ka.interpolate(mach, alt);
                Wi = obj.gainsThrustOff.Wi.interpolate(mach, alt);
                Kr = obj.gainsThrustOff.Kr.interpolate(mach, alt);
                AccMax = obj.gainsThrustOff.AccelMax.interpolate(mach, alt);
                rKx = obj.gainsThrustOff.rKx.interpolate(mach, alt);
                rKp = obj.gainsThrustOff.rKp.interpolate(mach, alt);
            end
            obj.Kdc = Kdc;
            obj.Ka = Ka;
            obj.Wi = Wi;
            obj.Kr = Kr;
            obj.rKx = rKx;
            obj.rKp = rKp;
        end

        function dpCmd = rollController(obj, phi, p)
            phiCmd = 0;
            dpCmd = 2.5 * obj.rKx * ( 5 * obj.rKp * ( phiCmd - phi) - p ); %[rad]
            dpCmd = obj.sat(dpCmd); % [rad]
        end

        function deltaCmd = pitchController(obj, errorPitch, q)
            deltaCmd = obj.Kr * (errorPitch - q);             
            deltaCmd = obj.sat(deltaCmd); % [rad]                        
        end

        function deltaCmd = yawController(obj, errorYaw, r)
            deltaCmd = obj.Kr * (errorYaw - r);  
            deltaCmd = obj.sat(deltaCmd); % [rad]
        end
        
        function errorDot = pitchDyn(obj, az, q, accCmd)
            az = az/obj.gravConst; % Convert Unit from m/s^2 -> gravity
            errorDot = obj.Wi * ( -obj.Ka * ( obj.Kdc * accCmd  - az ) - q );
        end
          
        function errorDot = yawDyn(obj, ay, r, accCmd)
            ay = ay/obj.gravConst; % Convert Unit from m/s^2 -> gravity
            errorDot = obj.Wi * ( obj.Ka * ( obj.Kdc * accCmd  - ay ) - r );
        end
          
        function [Kdc, Ka, omega_i, Kr] = polePlacing(obj, tranFcnCoff)
            a1 = tranFcnCoff.a1;
            a0 = tranFcnCoff.a0;
            b1 = tranFcnCoff.b1;
            b0 = tranFcnCoff.b0;
            c2 = tranFcnCoff.c2;
            c1 = tranFcnCoff.c1;
            c0 = tranFcnCoff.c0;
            g = obj.gravConst;
            Vm = obj.vehicle.speed;
            
            B = [obj.pole + 2 * obj.zeta * obj.omega - a1;...
                2 * obj.zeta * obj.omega * obj.pole + obj.omega^2 - a0;...
                obj.pole * obj.omega^2];
            A = [b1 0 -c2;...
                b0 b1 -c1;...
                0 b0 -c0];
            N = A\B;
            
            Kr = N(1);
            omega_i = N(2)/N(1);
            Ka = N(3)/N(2);
            Kdc = 1 + g/ ( Ka * Vm);
        end
        
        function out = commandShaping(obj, accCmd, accCdot, dt)
            accS = obj.shapedCmd;
            accSdot = obj.shapedCmdDot;
            accSddot = ( -0.0363 * accCdot + accCmd - 0.33 * accSdot - accS)/0.009;
            accSdot = accSddot * dt + obj.shapedCmdDot ;
            out = accS + obj.shapedCmdDot * dt;
            obj.shapedCmdDot = accSdot;
            obj.shapedCmd = out;
        end         

        function [scaledAcc1, scaledAcc2] = commandScaling(obj, threshold ,acc)
            accMag = vecnorm(acc, 2);
            threshold = min(threshold, accMag);
            accScaled = acc/( accMag +eps)* threshold;
            scaledAcc1 = accScaled(1);
            scaledAcc2 = accScaled(2);
        end
        
        function out = sat(obj, in)
            out = min(max(in, -obj.deltaCmdLim), obj.deltaCmdLim);
        end
    end
end
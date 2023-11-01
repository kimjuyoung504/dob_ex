classdef EarthAirEnv < handle
    properties
    end
    properties (Constant)
        airPressure = 101325  % atomospheric pressure [Pa]
        gasConstant = 287.053 % ideal gas constant [J/mol/K]
        gamma = 1.4           % Specific heat ratio
        G = 6.673e-11         % universal gravitational constant [Nm^2/kg^2]
        earthMass = 5.973e24  % Mass of the Earth   [kg]
        Re = 6370987.308      % Mean Earth Radius [m]
        T0 = 288.16           % Temp. at Sea Level [K]
        L = 0.0065            % Lapse Rate [K/m]
    end
    methods
        function obj = EarthAirEnv()
        end
        function T = Temperature(obj, alt)
            assert(alt >0, "Altitude must be larger than zero.");
            if alt < 11000
                T = obj.T0 - obj.L * alt;
            else
                T = 216;
            end
        end
        function P = Pressure(obj, alt)
            T = obj.Temperature(alt);
            if alt < 11000
                P = obj.airPressure * (T / 288.15 )^5.2559;
            else
                P = 22630 * exp( - 0.00015769 * ( alt - 11000 ) ) ;      
            end
        end
        function rho = density(obj, alt)
            T = obj.Temperature(alt);
            P = obj.Pressure(alt);
            rho = P / ( obj.gasConstant * T);
        end
        function a = soundSpeed(obj, alt)
            T = obj.Temperature(alt);
            a = sqrt( obj.gamma * obj.gasConstant * T);
        end
        function g = grav(obj)
            g = 0;
%             alt = -obj.vehicle.pos(3);
%             g = obj.G * obj.earthMass / (obj.Re + alt)^2;
        end
    end
end
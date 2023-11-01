% [Ref.1] Mauro Caresta, Vibration of a Free-Free beam.
classdef AeroelasticEffects < DynSystems
    properties (Hidden)
        modes % vector of the number of modes [-]
        E % Modulus of Elasticity [N/m^2]
        I % Moment of inertia [kg*m^2]
        L % Length of a beam [m]
        rho % Density of a beam [kg/m^3]
        area % Cross-section area [m^2]
        xi % Damping ratio [-]
        omega % Natural Frequency [1/rad]
        beta
        alpha
        mag % Initial states
        n_modes
        ignore = false
        dataNames
    end

    properties (Constant)
        beta_L = [4.7300408
                  7.8532046
                  10.995607839
                  14.1371654913
                  17.278759657399480]
    end

    methods
        function obj = AeroelasticEffects(name, modes, logOn)
            n_mode = length(modes);
            obj = obj@DynSystems(name, zeros(2*n_mode,1), logOn);
            obj.modes = modes;
            obj.n_modes = length(modes);
        end
        
        function obj = setParams(obj, E, I, L, rho, area, xi, mag, ignore)
            obj.ignore = ignore;
            obj.E = E;
            obj.I = I;
            obj.L = L;
            obj.rho = rho;
            obj.area = area;
            obj.beta = obj.beta_L(obj.modes)/L;
            obj.xi = xi;
            omega = obj.beta.^2 * ( E * I / (rho * area)).^(1/2);
            obj.mag = mag;
            for i = 1:length(obj.modes)
                q(i) = "q"+num2str(obj.modes(i));
                dqdt(i) ="dqdt "+num2str(obj.modes(i));
            end
            obj.omega = omega;
            obj.alpha = ( sin( obj.beta_L ) - sinh( obj.beta_L ) )...
                ./ ( cosh( obj.beta_L ) - cos( obj.beta_L ) );
            obj.dataNames = {q, dqdt, "zeta", "d\zeta /dt", "d\zeta dl"};
        end        
        
        function dsdt = dynEqns(obj, t, s, u)
            if obj.ignore == true
                u = zeros(obj.n_modes,1);
            end
            n_mode = length(obj.modes);
            B = [0; 1];
            dsdt = nan(2*n_mode,1);
            for i = 1:n_mode
                expA_t = obj.expAt(-t, obj.xi, obj.omega(i));
                dsdt(2*i-1:2*i) = expA_t * B * u(i);
            end
            [qs, qs_dot] = obj.get_q(t, s);

            if obj.logData
                obj.data.q = qs;
                obj.data.dqdt = qs_dot;
                obj.data.zeta = obj.elasticDisplacement(obj.L, qs);
                obj.data.dzetadt = obj.dzeta_dt(obj.L, qs_dot);
                obj.data.dzetadl = obj.dzeta_dl(obj.L, qs);
            end
        end

        function out = elasticDisplacement(obj, l, qs)
            phis = obj.phi(l);
            out = dot(phis, qs); % Eq. (2.8) in [Ref. 1]
        end
        
        function out = dzeta_dl(obj, l, qs)
            dphi_dl = obj.get_dphi_dl(l);
            out = dot(dphi_dl, qs);
        end

        function out =  dzeta_dt(obj, l, qs_dot)
            phis = obj.phi(l);
            out = dot(phis, qs_dot);
        end

        function out = phi(obj, l)
            l = l(:).';
            z = obj.alpha;
            out = z .* ( cosh( obj.beta * l ) + cos( obj.beta * l ) ) + sinh( obj.beta * l ) + sin( obj.beta * l );
        end

        function out = expAt(~, t, xi, omega)
            omega_d = omega * sqrt( 1 - xi^2);
            zeta = xi * omega;
            out = exp(-zeta * t)/omega_d *...
                   [omega_d * cos( omega_d * t ) + zeta * sin(omega_d * t), sin( omega_d * t );
                    -omega^2 * sin( omega_d * t ), omega_d * cos( omega_d * t) - zeta * sin( omega_d * t)];
        end

        function [qs, qs_dot] = get_q(obj, t, state)
            qs = nan(length(obj.omega), 1);
            qs_dot = nan(length(obj.omega), 1);
            for i = 1:length(obj.omega)
                x0 = [obj.mag;0];
                expAt = obj.expAt(t, obj.xi, obj.omega(i));
                temp = expAt*x0 + expAt * state(2*i-1:2*i);
                qs(i) = temp(1);
                qs_dot(i) = temp(2);
            end
        end

        function out = int_phi(obj, l)
            l = l(:).';
            z = obj.alpha;
            t = obj.beta * l;
            out = 1./obj.beta .* ( ( cosh( t ) - cos( t ) ) ...
                + z .* ( sinh( t ) + sin( t ) ) );
        end

        function out = int_phiSquare(obj, l)
            l = l(:).';
            z =  obj.alpha;
            t = obj.beta * l;
            b = obj.beta;
            out = l .* z.^2 ...
                  + z.^2 .* sin( 2 * t ) ./ (4 * b) ...
                  + z.^2 .* sinh( 2 * t ) ./ (4 * b) ...
                  + z.^2 .* cos( t ) .* sinh( t ) ./ b ...
                  + z.^2 .* sin( t ) .* cosh( t ) ./ b ...
                  - z .* cos(  2 * t ) ./ (2 * b) ...
                  + z .* cosh( 2 * t ) ./ (2 * b)...
                  + 2 * z .* sin( t ) .* sinh( t ) ./ b ...
                  - sin( 2 * t ) ./ (4 * b)...
                  + sinh( 2 * t ) ./ (4 * b) ...
                  - cos( t ) .* sinh( t ) ./ b...
                  + sin( t ) .* cosh( t ) ./ b;
        end

        function out = get_dphi_dl(obj, l)
            b = obj.beta;
            z = obj.alpha;
            out = b .* ( cosh(b .* l) + cos(b .* l)  + z .* ( sinh( b .* l) - sin( b .* l ) ) );
        end

        function out = generalizedForce(obj, force)
            out = force/obj.L * obj.phi(obj.L);
        end

        function out = generalizedMass(obj, mass)
            out = mass/obj.L * obj.int_phiSquare(obj.L);
        end

        function alpha = inducedIncidesAngles(obj, u, alpha, l, q, q_dot)
            alpha = alpha - obj.dzeta_dl(l, q) - obj.dzeta_dt(l, q_dot)/u;
        end
    end

    methods (Static)
        function test()
%             set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.5,'defaultfigurewindowstyle','docked');
%             set(0,'defaultfigurecolor','white');
%             set(0, 'DefaultLineLineWidth', 1);
            % Data of a beam
            E = 2.1*1e11; % 2.1 x 10^11 N/m
            L = 1.275; % m
            rho = 7800; % Kg m^-3
            b = 0.075; % m
            h = 0.01; % m
            area = h*b; % m^2
            I = b*h^3/12;
            xi = 0.0001;

            % Set system
            mag = 0.01;
            modes = [1,2,3,4,5];
            n_mode = length(modes);
            logOn = true;
            sys = AeroelasticEffects("test", modes, logOn).setParams(E, I, L, rho, area, xi, mag);
            
            % Check mode shape
            l = 0:0.01:L;
            phis = sys.phi(l);
            int_phis = sys.int_phi(l);
            int_phiSq = sys.int_phiSquare(l);

            figure(1); hold on;
            for i = 1:length(modes)
                plot(l, phis(modes(i), :))
            end
            legend("1st mode", "2nd mode", "3rd mode", "4th mode","5th mode")

            % Solve the dynamical system
            input = @(t) zeros(5,1);
            t0 = 0;
            dt = 0.0025;
            tf = 2.0;
            tspan = t0:dt:tf;
            sim = Simulator(sys).propagate(tspan, input);
            sim.report();
            log = sim.log;

            % Plot results
            tspan = 0:dt:tf;
            state_mode = nan(length(tspan), length(modes));
            for i = 1:length(tspan)
                for j = 1:length(modes)
                    expAt = sys.expAt(tspan(i), xi, sys.omega(j));
                    s = expAt * [mag; 0];
                    state_mode(i, j) = s(1);
                end
            end
            log.q.subplots(2);
            log.zeta.plot(3);
            figure(2); hold on;
            for i = 1:length(modes)
                subplot(length(modes),1,i)
                plot(tspan, state_mode(:, i),"--")
                legend("Numerical", "Analytic")
            end
            figure(4); hold on;
            for i = 1:n_mode
                plot(l, int_phis(i, :))
            end
            legend("1st mode", "2nd mode", "3rd mode", "4th mode","5th mode")
            
            figure(5); hold on;
            for i = 1:n_mode
                plot(l, int_phiSq(i, :))
            end            

            log.dqdt.subplots(6);
            log.zeta.plot(7);
            log.dzetadt.plot(8);
            log.dzetadl.plot(9);
        end
    end
end
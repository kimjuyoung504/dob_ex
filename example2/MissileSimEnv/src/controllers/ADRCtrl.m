classdef ADRCtrl < DynSystems
    properties
        dyn_order
        omega_o
        omega_c
        b_hat
        u_lim
        T_d
    end
    methods
        function obj = ADRCtrl(name, dyn_order, logOn)
            x0 = zeros(dyn_order+1,1);
            obj = obj@DynSystems(name, x0, logOn);
            obj.dyn_order = dyn_order;
        end
        
        function obj = setParams(obj, omega_o, omega_c, b_hat, u_lim, T_d)
            obj.omega_o = omega_o;
            obj.omega_c = omega_c;
            obj.b_hat = b_hat;
            obj.u_lim = u_lim;
            obj.T_d = T_d;
        end
        function dsdt = dynEqns(obj, ~, s, y, u)
            A = obj.get_A();
            B = obj.get_B();
            C = obj.get_C();            
            ko = obj.get_Ko();

            x_hat_dot =  A * s + B * u +...
                         ko * ( y - dot(C, s));
            dsdt = x_hat_dot;
        end

        function u = get_command(obj, s_hat, t, x_d)
            x_hat = s_hat(1:end-1);
            d_hat = s_hat(end);
            Kc = obj.get_Kc();

            u = 1/obj.b_hat * ( Kc(1) * x_d - dot(Kc, x_hat) - d_hat);
            u = min(max(u, obj.u_lim(1)), obj.u_lim(2)); % saturation
            if t  <  obj. T_d % Anti-peaking
                u = 0;
            end
        end

        function A = get_A(obj)
            do = obj.dyn_order;
            temp = [eye(do); zeros(1, do)];
            A = [zeros(do+1,1), temp];
        end

        function B = get_B(obj)
            do = obj.dyn_order;
            B = [zeros(do-1, 1); obj.b_hat; 1];
        end

        function C = get_C(obj)
            do = obj.dyn_order;
            C = [1, zeros(1, do)];
        end

        function Kc = get_Kc(obj)
            do = obj.dyn_order;
            oc = obj.omega_c;
            temp = poly(-oc * ones(do, 1));
            Kc = flip(temp(2:end));
        end
        
        function Ko = get_Ko(obj)
            do = obj.dyn_order;
            temp = poly(-obj.omega_o * ones(do+1,1) );
            Ko = temp(2:end).';
        end
    end

    methods (Static)
        function test()
            f_star = @(t,x) -4.0*x(4) -6.0*x(3) -4.0*x(2) -x(1);
            b_star = 1.0;
            d_star = @(t) 70.0*heaviside(t-10);

            n = 4;
            b_hat = 0.8;
            x_d = 1;

            omega_c = 5;
            omega_o = 50;
            u_min = -100;
            u_max = 100;
            std = 1e-17;
            u_lim = [u_min, u_max];
            Td = 0.1;
            logOn = true;

            adrc = ADRCtrl("adrc", n, logOn).setParams(omega_o, omega_c, ...
                b_hat, u_lim, Td);
            
            tf = 20;
            s0 = zeros(9, 1);
            tspan = [0, tf];
            
            [t, s]= ode45(@f, tspan, s0);
            
            plot(t, s(:, 1))
            function dsdt = f(t, s)
                x = s(1:4);
                y = x(1) + std*randn();
                x_hat = s(5:end);
                u = adrc.get_command(x_hat, t, x_d);
                z_hat_dot = adrc.dynEqns(t, x_hat, y, u);
                x4_dot = f_star(t, x) + b_star *  u + d_star(t);
                x_dot = [x(2:4); x4_dot];
                dsdt = [x_dot; z_hat_dot];
            end
        end
    end
end
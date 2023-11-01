%Gernot Herbst, "A Simulative Study on Active Disturbance Rejection Control
%(ADRC) as a Control Tool for Practitioners," arXiv2019, doi:arXiv:1908.04596.
classdef ADRCIRL < DynSystems
    properties
        oc %omega_c = bandwidth
        op_g
        b0
        order
        r
        limit
        Td 
        Q
        R
        u0
    end
    methods 
        function obj = ADRCIRL(name, x0)
            obj = obj@DynSystems(name, x0, true);
        end

        function obj = setParams(obj, b0, oc, order, op_g, r, limit, Td, Q, R)
            obj.b0 = b0;
            obj.oc = oc;
            obj.order = order;
            obj.r = r;
            obj.op_g = op_g;
            obj.limit = limit;
            obj.Td = Td;
            obj.R = R;
            obj.Q = Q;
        end

        function out = get_A(obj)
            out = [zeros(obj.order, 1), eye(obj.order); zeros(1, obj.order + 1)];
        end

        function out = get_B(obj)
            out = [zeros(obj.order - 1, 1); obj.b0; 0];
        end
        function out = get_C(obj)
            out = [1, zeros(1, obj.order)];
        end

        function out = get_L(obj)
            oc_op = obj.op_g * obj.oc;
            param = poly(-oc_op * ones(1, obj.order + 1));
            out = param(2:end).';
        end
        
        function out = Value(obj)
        end

        function u = controller(obj, t, x, u0)
            u = 1/obj.b0 * (u0 - x(end));
            u = min(max(u, -obj.limit), obj.limit);
            if t < obj.Td
                u = 0;
            end
        end
       

        function dx = dynEqns(obj, ~, x, input)
            u = input(1:end-1);
            y = input(end);
            A = obj.get_A();
            L = obj.get_L();
            B = obj.get_B();
            C = obj.get_C();

            dx = (A -  L * C) * x + B * u + L * y;
        end 
    end

    methods (Static)
        function test()
            f_star = @(t,x) -4.0*x(4) -6.0*x(3) -4.0*x(2) -x(1);
            b_star = 1.0;
            d_star = @(t) 70.0*heaviside(t-10);

            order = 4;
            b0 = 0.8;
            r = 1;

            oc = 5;
            op_g = 10;
            std = 1e-17;
            limit = 100;
            Td = 0.1;
            x0 = zeros(order+1, 1);
            adrc = ADRC2("adrc", x0).setParams(b0, oc, order, op_g, r, limit, Td);

            tf = 20;
            s0 = [zeros(order, 1); x0];
            tspan = [0, tf];

            [t, s]= ode45(@f, tspan, s0);

            plot(t, s(:, 1)); hold on;
            plot(t, r*ones(length(t), 1))
            legend("Actual", "Reference")
            function dsdt = f(t, s)
                x = s(1:4);
                y = x(1) + std*randn();
                x_hat = s(5:end);
                u = adrc.controller(t, x_hat);
                z_hat_dot = adrc.dynEqns(t, x_hat, [u, y]);
                x4_dot = f_star(t, x) + b_star *  u + d_star(t);
                x_dot = [x(2:4); x4_dot];
                dsdt = [x_dot; z_hat_dot];
            end
        end
    end

end
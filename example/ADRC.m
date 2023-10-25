classdef ADRC < handle
    properties
        oc = 1; %omega_c = bandwidth
        op_g = 10;
        b0 = 1;
        order = 2;
        r = 1;
        limit = 20;
    end
    methods
        function obj = ADRC()
        end
        function set_parameter(obj, b0, oc, order, op_g, r)
            obj.b0 = b0;
            obj.oc = oc;
            obj.order = order;
            obj.r = r;
            obj.op_g = op_g;
        end
        function out = get_A(obj)
            out = [zeros(obj.order, 1), eye(obj.order); 0, zeros(1, obj.order)];
        end
        function out = get_B(obj)
            out = [zeros(obj.order - 1, 1); obj.b0; 0];
        end
        function out = get_C(obj)
            out = [1, zeros(1, obj.order)];
        end
        function out = get_L(obj)
            oc_op = obj.op_g * obj.oc;
            root = oc_op * ones(1, obj.order + 1);
            out = flip(poly(root));
        end
        function out = K(obj)
            root = obj.oc * ones(1, obj.order + 1);
            out = flip(polu(root));
        end
        
        function u = controller(obj, x)
            u = 1/obj.b0 * (K(1) * obj.r + dot(K(obj) * x(1:obj.order)) - x(obj.order+1));
        end
        function u_l = u_limit(obj, x)
            u_l = min(max(controller(obj, x), -deg2rad(obj.limit)), deg2rad(obj.limit));
        end
        function dx = esti_dyn(obj, x, y)
            A = obj.get_A();
            L = obj.get_L();
            B = obj.get_B();
            C = obj.get_C();
            input = [controller(obj, x), y];
            dx = (A -  L * C) * x + [B, L] * input;
        end 
    end

end

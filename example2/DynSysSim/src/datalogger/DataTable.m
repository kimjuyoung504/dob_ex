classdef DataTable < handle
    properties
        dim
        points
        data
        Fcn
    end
    properties (Hidden)
        method = 'linear'     
        grid
    end  
    methods
        function obj = DataTable(data, varargin)
           obj.dim = numel(varargin);           
           obj.points = varargin;
           obj.data = data;
           obj.grid = cell(obj.dim, 1);     
           if obj.dim == 1
               obj.Fcn = griddedInterpolant(varargin{:}, data, obj.method);
           else
               [obj.grid{:}] = ndgrid(varargin{:});   
               obj.Fcn = griddedInterpolant(obj.grid{:}, data);
           end
        end
        
        function out = interpolate(obj, varargin)
            if numel(varargin) ~= obj.dim
                error('Dimension of quary point and table does not match');
            end
            out = obj.Fcn(varargin{:});
        end              
    end
    
    methods (Static)
        function test()
            disp('Check if the interpolate method is well operated.');
            %% 1-d example
            x = 1:10;
            V = rand(1, 10);
            self = DataTable(V, x);
            out = self.interpolate([5 5.5 7]);
            disp(['1-D case : ',num2str(out)])
            %% 2-d example
            % create example data 
            x = 1:10;
            y = 1:5;
            V = rand(10, 5);
            
            self = DataTable(V, x, y);
            out = self.interpolate([1, 2, 3],[1, 4, 5]);
            disp(['2-D case : ',num2str(out)])
            
            %% 3-d example
            z = x;
            V = rand(10, 5, 10);
            
            self = DataTable(V, x, y, z);
            out = self.interpolate([1, 2, 4],[1, 3, 4],[1, 4, 4]);
            disp(['3-D case : ',num2str(out)])
        end
    end
end
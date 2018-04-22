classdef Rozenbrock
    methods
        function obj = Rozenbrock()
        end
        
        function z = Func(obj,x)
            z = 100*(x(2,:) - x(1,:).^2).^2 + (1-x(1,:)).^2; 
        end
        
        function z = FuncMesh(obj,x,y)
            z = 100*(y - x.^2).^2 + (1-x).^2; 
        end
        
        function z = Deriv(obj,x)
           z = ([-400.*(x(2,:) - x(1,:).^2).*x(1,:) - 2.*(1-x(1,:)) ; 
               200*(x(2,:) - x(1,:).^2)])';

        end
    end
    
end


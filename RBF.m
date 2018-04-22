classdef RBF < RBFbase
    methods
        function obj = RBF(x, y, varargin)
            obj@RBFbase(x, y, varargin);
        end
    end
    
    methods (Access = protected)
        function [obj, baseFuncName] = ParseVarargin(obj,varargin)
            [XDim, XCount] = size(obj.nodes);
            obj.constant = (prod(max(obj.nodes') - min(obj.nodes'))/XCount)^(1/XDim);
            obj.smooth   = 0;
            varargin = varargin{1,1};
             varargin = varargin{1};
            if length(varargin) > 0
                baseFuncName = varargin{1};
            else
                baseFuncName = 'linear'; 
            end
        end
        
        function [obj, A] = Assemble(obj)
            [dim, n] = size(obj.nodes);
            A = zeros(n, n);
            for i = 1:n
                for j = 1:i
                    distant = norm(obj.nodes(:,i) - obj.nodes(:,j));
                    temp = feval(obj.baseFunction.Func, distant, obj.constant);
                    A(i,j) = temp;
                    A(j,i) = temp;
                end
                A(i,i) = A(i,i) - obj.smooth;
            end
        end  
        
        
         function y = FormY(obj,y)
         end

    end
end


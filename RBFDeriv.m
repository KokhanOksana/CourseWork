classdef RBFDeriv < RBFbase   
    properties (Access = protected)
        derivValues,
        aditionalNodes
    end
    methods
        function obj = RBFDeriv(x, y, varargin)
            obj@RBFbase(x, y, varargin);
        end
    end
    
   methods (Access = protected)
        function [obj, baseFuncName] = ParseVarargin(obj,varargin)
            [XDim, XCount] = size(obj.nodes);
            obj.constant = 2;%(prod(max(obj.nodes') - min(obj.nodes'))/XCount)^(1/XDim);
            obj.smooth   = 0;
            varargin = varargin{1,1};
             varargin = varargin{1};
            if length(varargin) > 0
                baseFuncName = varargin{1};
            else
                baseFuncName = 'linear'; 
            end
            if length(varargin) > 1
                obj.derivValues = varargin{2};
            end
             if length(varargin) > 2
                obj.aditionalNodes = varargin{3};
            end
        end
        
        function [obj,A] = Assemble(obj)
           [dim, n] = size(obj.nodes);
           obj.nodes = [obj.nodes obj.aditionalNodes];
            A = zeros(n*(dim+1),n + size(obj.aditionalNodes,2));
            for i = 1:n
                for j = 1:size(obj.nodes,2)
                    distant = norm(obj.nodes(:,i) - obj.nodes(:,j));
                    A(i,j) = feval(obj.baseFunction.Func, distant, obj.constant);
                    deriv =  obj.baseFunction.normDeriv(obj.nodes(:,i), obj.nodes(:,j));
                    for k = 1: dim
                        A(i+k*n,j) = feval(obj.baseFunction.Deriv,deriv(k),obj.constant);
                    end
                end
                for k = 1: dim
                    A(i+k*n,i) = A(i,i) - obj.smooth;
                end
            end
        end
        
         function y = FormY(obj,y)
            dim = size(obj.derivValues,2);
            for i = 1:dim;
                y = [y ;obj.derivValues(:,i)];  
            end
         end
    end
end

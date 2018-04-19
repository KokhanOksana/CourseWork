classdef RBF
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        smooth 
        constant 
        nodes
        derivValues
        baseFunction
        weightCoeff
        aditionalNodes
    end
    
    %varargin(1) - BaseFunctionName
    %varargin(2)- DerivValues
    methods (Access = public)
        function obj = RBF(x, y, varargin)
            obj.nodes = x;
            
            [obj,baseFuncName] = obj.ParseVarargin(varargin);
            obj.baseFunction = BaseFunction(baseFuncName);
            
            if length(varargin) > 1
                [obj,A] = obj.AssembleGrad();
                y = obj.AddDerivatives(y);
            else
                A = obj.Assemble();
            end
            obj.weightCoeff =  A\y;
        end
        
        function [f] = Interpolate(obj, x)
            n = size(obj.nodes,2);
            nPoints = size(x,2);

            f = zeros(1, nPoints);
            r = zeros(1, n);

            for i = 1:1:nPoints
                r =  (x(:,i)*ones(1,n)) - obj.nodes;
                r = sqrt(sum(r.*r, 1));
                f(i) = sum(obj.weightCoeff'.*feval(obj.baseFunction.Func, r, obj.constant));
            end
        end
    end 
    methods (Access = private)
        function [obj, baseFuncName] = ParseVarargin(obj,varargin)
            [XDim, XCount] = size(obj.nodes);
            obj.constant = (prod(max(obj.nodes') - min(obj.nodes'))/XCount)^(1/XDim);
            obj.smooth   = 0;
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
        
        function [A] = Assemble(obj)
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
        function [obj,A] = AssembleGrad(obj)
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
        
         function y = AddDerivatives(obj,y)
            dim = size(obj.derivValues,2);
            for i = 1:dim;
                y = [y ;obj.derivValues(:,i)];  
            end
         end

    end
end


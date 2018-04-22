classdef (Abstract) RBFbase
    
   properties (Access = protected)
        smooth 
        constant 
        nodes
        baseFunction
        weightCoeff
   end
    
    methods
     function obj = RBFbase(x, y, varargin)
            obj.nodes = x;
            
            [obj,baseFuncName] = obj.ParseVarargin(varargin);
            obj.baseFunction = BaseFunction(baseFuncName);
            [obj,A] = obj.Assemble();
            y = obj.FormY(y);
           
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
                w = obj.weightCoeff';
                p = feval(obj.baseFunction.Func, r, obj.constant);
                f(i) = sum(w .* p);
            end
     end
    end
     methods (Access = protected,Abstract)
      [obj, baseFuncName] = ParseVarargin(obj,varargin)
      [obj,A] = Assemble(obj)
      y = FormY(obj,y)   
     end
end

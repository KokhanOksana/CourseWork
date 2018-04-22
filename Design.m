classdef Design
    properties
        x, f, fDeriv, aditionalNodes, 
        count, aditionalCount
    end
    
    methods
        function obj = Design()
            obj.count = 0;
        end
        
        function obj = LH(obj, count, bound, func)
            global callCount; 
            callCount = callCount + count;
            obj.count = count;
            obj.x = obj.LHPoints(count, bound);
            obj.f = func.Func(obj.x');
            obj.fDeriv = func.Deriv(obj.x');
        end
        
        function x = LHPoints(obj, count, bound)
            global dim;
            xLH = lhsdesign(count,dim);
            x = zeros(count, dim);
            for j = 1: count
                for p = 1:dim
                    d = bound.b(p) - bound.a(p);
                    x(j,p) = bound.a(p) + d * xLH(j,p);
                end
            end
        end
        
        function obj = New(obj, x, f, fDeriv)
            obj.count = size(x,1);
            obj.x = x;
            obj.f = f;
            obj.fDeriv = fDeriv;
        end
        
        function obj = Copy(obj, design)
            obj.count = design.count;
            obj.x = design.x;
            obj.f = design.f;
            obj.fDeriv = design.fDeriv;
        end
        
         function [obj] = AddFromDesign(obj, design, designIndex)
            obj.count = obj.count+1;
            obj.x(obj.count,:) = design.x(designIndex,:);
            obj.f(obj.count) = design.f(designIndex);
            obj.fDeriv(obj.count,:) = design.fDeriv(designIndex,:);
         end
       
    end
    
end


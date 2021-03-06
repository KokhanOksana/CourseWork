classdef MetaModelDeriv < MetaModel 
    methods
        function obj = MetaModelDeriv(fun, eps, nodeCount)
            obj@MetaModel(fun, eps, nodeCount);
        end
    end
    
    methods(Access = protected)
        function [ xMin, fMin] = SubRegionMin(obj)
            aditionalNodes = obj.design.LHPoints(obj.nodeCount/2, obj.bound);
            rbf = RBFDeriv(obj.design.x', obj.design.f', 'multiquadric', obj.design.fDeriv', aditionalNodes'); 
            min_opt = optimset('Display','off');
            [xMin, fMin] = fmincon(@rbf.Interpolate, obj.xMin,[],[],[],[], obj.bound.a,obj.bound.b,[], min_opt);
        end  
        
        function flag = BreakCondition(obj,radius, xMinPrev, fMinPrev)
            global dim
            flag = 0;
             if(obj.iteration > 100)
                 display('step abuse');
                 flag = 1;
             end;
             if(max(abs( xMinPrev - obj.xMin))< obj.eps)
                 display('points closeness');
                 flag = 1;
             end
             xMinDeriv = obj.func.Deriv(obj.xMin)
             if(sum(abs(xMinDeriv)) < obj.eps)
                 display('deriv = 0');
                 flag = 1;
             end
        end
        
        function deriv =  CloseDesignPointDeriv(obj, x)
            d = inf;
            for i = 1: size(obj.design.x,1)
                if(distant(obj.design.x(i,:),x) < d)
                    d = distant(obj.design.x(i,:),x);
                    deriv = obj.design.fDeriv(i,:);
                end
            end
            deriv = sum(deriv);
        end
    end
end


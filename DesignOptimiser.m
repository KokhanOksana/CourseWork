classdef DesignOptimiser
    properties
        bound, prevBound,
        design, prevDesign,
    end
    
    methods
        function obj = DesignOptimiser(bound, prevBound, design, prevDesign)
            obj.bound = bound;
            obj.prevBound = prevBound;
            obj.design = design;
            obj.prevDesign = prevDesign;
        end
        
        function [newDesign] = Optimise(obj,nodeCount, fun)
            newDesign = Design();
            for k = 1:size(obj.prevDesign.x,1)
                [takePrev, takeDesign] = obj.Take(k);
                if(takePrev == 1)
                    newDesign = newDesign.AddFromDesign(obj.prevDesign, k);
                end
                if(takeDesign == 1)
                    newDesign = newDesign.AddFromDesign(obj.design, k);
                end
            end;
            if(newDesign.count < nodeCount)
                aditional = Design();
                aditional = aditional.LH(nodeCount - newDesign.count, obj.bound, fun);
                for k = 1:size(aditional.x)
                   newDesign = newDesign.AddFromDesign(aditional, k); 
                end
            end
        end
        
        function [takePrev, takeDesign] = Take(obj, l)
            global dim
            takePrev = 1;
            takeDesign = 0;
            for k = 1:dim
                if((obj.prevDesign.x(l,k) < max(obj.bound.a(k),obj.prevBound.a(k))) || (obj.prevDesign.x(l,k) > min(obj.bound.b(k),obj.prevBound.b(k))))
                    takePrev = 0;
                end
                if(l <= size(obj.design.x,1)) 
                    if((obj.design.x(l,k) > obj.prevBound.b(k)) || (obj.design.x(l,k) < obj.prevBound.a(k)))
                        takeDesign = 1;
                    end
                end
            end
        end
    end
    
end


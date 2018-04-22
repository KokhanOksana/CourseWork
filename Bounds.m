classdef Bounds
    properties
        a, b, sizes,
        aStart, bStart,
        minRadius
    end
    
    methods
        function obj = Bounds(a, b)
            obj.a = a;
            obj.b = b;
            obj.sizes = b - a;
            obj.aStart = a;
            obj.bStart = b;
            obj.minRadius = max(b - a)/10;
        end
        
        function  obj = Bound(bound)
            obj.a = bound.a;
            obj.b = bound.b;
            obj.sizes = bound.sizes;
            obj.aStart = bound.aStart;
            obj.bStart = bound.bStart;
        end
        
        function obj = New(obj, xMin, radius)
        global dim;
            obj.a = xMin' - radius*ones(1,dim);
            obj.b = xMin' + radius*ones(1,dim);
            obj = obj.Cut();

            plot(xMin(1), xMin(2), '*r');hold on;
            obj.PlotBound();
        end
        
        function [] = PlotBound(obj)
            global dim;            
            linex = obj.a(1):0.01:obj.b(1);
            lineLengthx = size(linex,2);
            liney = obj.a(2):0.01:obj.b(2);
            lineLengthy = size(liney,2);
            
            plot( linex, obj.b(2)*ones(1,lineLengthx) ,'-g');hold on;
            plot( linex, obj.a(2)*ones(1,lineLengthx) ,'-g');hold on;
            plot( obj.a(1)*ones(1,lineLengthy), liney ,'-g');hold on;
            plot( obj.b(1)*ones(1,lineLengthy), liney ,'-g');hold on;
        end
        
        function obj = Cut(obj)
            global dim;    
            for k = 1 :dim
                if(obj.a(k) < obj.aStart(k))
                    obj.a(k) = obj.aStart(k);
                end
                if(obj.b(k) > obj.bStart(k))
                    obj.b(k) = obj.bStart(k);
                end
            end;
        end
        
    end
    
end


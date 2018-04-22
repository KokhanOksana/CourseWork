classdef MetaModel
    properties
        iteration, xMin, fMin,
    end
    properties(Access = protected)
        func, eps, nodeCount,
        bound, design
    end
    
    methods
        function obj = MetaModel(fun, eps, nodeCount)
            obj.func = fun;
            obj.eps = eps;
            obj.nodeCount = nodeCount;
            obj.iteration = 0;
        end
        
        function obj = Min( obj, x0, bound0 )
            obj = obj.Init(x0, bound0);
            obj.design = Design();
            obj.design = obj.design.LH(obj.nodeCount, obj.bound, obj.func);
            dg = inf;
            df = inf;
            obj.ShowResult();
            while(1)
                obj.iteration = obj.iteration + 1;
                
                xMinPrev = obj.xMin;
                fMinPrev = obj.fMin;
                [obj.xMin, obj.fMin] = obj.SubRegionMin();

                radius = distant(obj.xMin, xMinPrev);
                obj = obj.NewBoundDesign(radius); 
                
                if(obj.BreakCondition(radius, xMinPrev, fMinPrev))
                    break
                end
                obj.ShowResult();
            end;
        end
        
       function obj = ShowResult(obj)
           global callCount;
           names = {'iteration','callCount', 'xMin1','xMin2','fMin'};
           display(table([obj.iteration; callCount;  obj.xMin(); obj.fMin], 'RowNames', names));
        end
        
    end
    
    methods(Access = protected)
         function obj = Init(obj, x0, bound0)
            obj.bound = bound0;
            obj.xMin = x0;
            obj.fMin = obj.func.Func(obj.xMin);
        end
        
        function [ xMin, fMin] = SubRegionMin(obj)
            rbf = RBF(obj.design.x', obj.design.f', 'gaussian'); 
            min_opt = optimset('Display','off');
            [xMin, fMin] = fmincon(@rbf.Interpolate, obj.xMin,[],[],[],[], obj.bound.a,obj.bound.b,[], min_opt);
        end
        
        function obj = NewBoundDesign(obj,radius)
            if(radius < obj.bound.minRadius)
                radius = obj.bound.minRadius;
            end
            %prevBound = Bound(obj.bound);
            obj.bound = obj.bound.New(obj.xMin, radius);
            %prevDesign = Design().Copy(obj.design);
            obj.design = obj.design.LH(obj.nodeCount, obj.bound, obj.func);
%             optimiser = DesignOptimiser(obj.bound, prevBound, obj.design, prevDesign);
%             obj.design = optimiser.Optimise(obj.nodeCount, obj.func);
        end
        
        function flag = BreakCondition(obj,radius, xMinPrev, fMinPrev)
            flag = 0;
             if(obj.iteration > 100)
                 display('step abuse');
                 flag = 1;
             end;
             if(max(abs( xMinPrev - obj.xMin))< obj.eps)
                 display('points closeness');
                 flag = 1;
             end
             if(abs((obj.fMin - fMinPrev) / radius) < obj.eps )
                  display('deriv = 0');
                 flag = 1;
             end
        end
        
    end
    
    end


classdef RBF
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        Smooth 
        Constant 
        Nodes
        DerivValues
        BaseFunctionName
        BaseFunction
        BaseFunctionDeriv
        WeightCoeff
    end
    
    methods (Access = public)
        function obj = RBF(x, y, varargin)
            obj.Nodes = x;
            [XDim, XCount] = size(obj.Nodes);
            
            obj = obj.ParseVarargin(varargin)
            obj = obj.SetBaseFunction();
            
            A = obj.Assemble();
            obj.WeightCoeff =  A\y';
        end
        
        function [f] = Interpolate(obj, x)
            [dim              n] = size(obj.Nodes);
            [dimPoints  nPoints] = size(x);

            f = zeros(1, nPoints);
            r = zeros(1, n);

            for i = 1:1:nPoints
                r =  (x(:,i)*ones(1,n)) - obj.Nodes;
                r = sqrt(sum(r.*r, 1));
                f(i) = sum(obj.WeightCoeff(1:n)'.*feval(obj.BaseFunction, r, obj.Constant));
            end
        end
    end 
    methods (Access = private)
        function obj = ParseVarargin(obj,varargin)
            [XDim, XCount] = size(obj.Nodes);
            if length(varargin) > 0
                obj.BaseFunctionName = varargin{1};
            else
                obj.BaseFunctionName = 'linear'; 
            end
            if length(varargin) > 1
                obj.Constant = varargin{2};
            else
                obj.Constant = (prod(max(obj.Nodes') - min(obj.Nodes'))/XCount)^(1/XDim);
            end
            if length(varargin) > 2
                obj.Smooth = varargin{3};
            else
                obj.Smooth   = 0;
            end
            if length(varargin) > 3
                obj.DerivValues = varargin{4};
            else
                obj.DerivValues   = [];
            end
        end
        function obj = SetBaseFunction(obj)
            switch lower(obj.BaseFunctionName{1})
                  case 'linear'          
                    obj.BaseFunction   = @obj.rbfphi_linear;
                  case 'cubic'
                    obj.BaseFunction   = @obj.rbfphi_cubic;
                  case 'multiquadric'
                    obj.BaseFunction   = @obj.rbfphi_multiquadrics;
                  case 'thinplate'
                    obj.BaseFunction   = @obj.rbfphi_thinplate;
                  case 'gaussian'
                    obj.BaseFunction   = @obj.rbfphi_gaussian;
                  case 'cubicspline'
                    obj.BaseFunction   = @obj.rbfphi_cubicspline; 
                otherwise
                    obj.BaseFunction   = @obj.rbfphi_linear;
            end
        end
        
        function [A] = Assemble(obj)
            [dim, n] = size(obj.Nodes);
            A = zeros(n, n);
            for i = 1:n
                for j = 1:i
                    distant = norm(obj.Nodes(:,i) - obj.Nodes(:,j));
                    temp = feval(obj.BaseFunction, distant, obj.Constant);
                    A(i,j) = temp;
                    A(j,i) = temp;
                end
                A(i,i) = A(i,i) - obj.Smooth;
            end
        end
        function [A] = AssembleGrad()
           [dim, n] = size(obj.Nodes);
            A = zeros(n, n*(dim+1));
            for i = 1:n
                for j = 1:i
                    distant = norm(obj.Nodes(:,i) - obj.Nodes(:,j));
                    temp = feval(obj.BaseFunction, distant, obj.Constant);
                    A(i,j) = temp;
                    A(j,i) = temp;
                end
                A(i,i) = A(i,i) - obj.Smooth;
            end
            for k = 1:dim
                
            end
        end
        
        %**************************************************************************
        % Radial Base Functions
        %************************************************************************** 
        function u=rbfphi_linear(obj, r, const)
            u=r; end
        function u=rbfphi_cubic(obj, r, const)
            u=r.*r.*r;end
        function u=rbfphi_gaussian(obj, r, const)
            u=exp(-0.5*r.*r/(const*const));end
        function u=rbfphi_multiquadrics(obj, r, const)
            u=sqrt(1+r.*r/(const*const));end
        function u=rbfphi_thinplate(obj, r, const)
            u=r.*r.*log(r+1);end
        function u=rbfphi_cubicspline(obj, r, const)
            u= (r.*r + const*const).^(3/2);end
        
        % Radial Base Functions Derivatives
        function u=rbfphi_deriv(obj, r, const)
            u=r; end
        function u=rbfphi_cubic_deriv(obj, r, const)
            u=r.*r.*r;end
        function u=rbfphi_gaussian_deriv(obj, r, const)
            u=exp(-0.5*r.*r/(const*const));end
        function u=rbfphi_multiquadrics_deriv(obj, r, const)
            u=sqrt(1+r.*r/(const*const));end
        function u=rbfphi_thinplate_deriv(obj, r, const)
            u=r.*r.*log(r+1);end
        function u=rbfphi_cubicspline_deriv(obj, r, const)
            u= (r.*r + const*const).^(3/2);end
    end
end


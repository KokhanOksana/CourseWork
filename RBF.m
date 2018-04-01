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
    
    %varargin(1) - BaseFunctionName
    %varargin(2) - Constant
    %varargin(3) - Smooth
    %varargin(4) - DerivValues
    methods (Access = public)
        function obj = RBF(x, y, varargin)
            obj.Nodes = x;
            [XDim, XCount] = size(obj.Nodes);
            
            obj = obj.ParseVarargin(varargin)
            obj = obj.SetBaseFunction();
            
            %A = obj.Assemble();
            
            A = obj.AssembleGrad();
            y = obj.AddDerivatives(y);
            
            obj.WeightCoeff =  A\y;
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
            varargin = varargin{1};
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
        
        function y = AddDerivatives(obj,y)
            [n, dim] = size(obj.DerivValues);
            for i = 1:dim;
                y = [y ;obj.DerivValues(:,i)];  
            end
        end
        
        function obj = SetBaseFunction(obj)
            switch lower(obj.BaseFunctionName)
                  case 'linear'          
                    obj.BaseFunction   = @obj.rbfphi_linear;
                    obj.BaseFunctionDeriv = @obj.rbfphi_linear_deriv;
                  case 'cubic'
                    obj.BaseFunction   = @obj.rbfphi_cubic;
                    obj.BaseFunctionDeriv =  @obj.rbfphi_cubic_deriv;
                  case 'multiquadric'
                    obj.BaseFunction   = @obj.rbfphi_multiquadrics;
                    obj.BaseFunctionDeriv = @obj.rbfphi_multiquadrics_deriv;
                  case 'thinplate'
                    obj.BaseFunction   = @obj.rbfphi_thinplate;
                    obj.BaseFunctionDeriv =  @obj.rbfphi_thinplate_deriv;
                  case 'gaussian'
                    obj.BaseFunction   = @obj.rbfphi_gaussian;
                    obj.BaseFunctionDeriv = @obj.rbfphi_gaussian_deriv;
                  case 'cubicspline'
                    obj.BaseFunction   = @obj.rbfphi_cubicspline;
                    obj.BaseFunctionDeriv = @obj.rbfphi_cubicspline_deriv;
                otherwise
                    obj.BaseFunction   = @obj.rbfphi_linear;
                    obj.BaseFunctionDeriv =  @obj.rbfphi_linear_deriv;
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
        
        function [A] = AssembleGrad(obj)
           [dim, n] = size(obj.Nodes);
            A = zeros(n*(dim+1),n);
            for i = 1:n
                for j = 1:i
                    distant = norm(obj.Nodes(:,i) - obj.Nodes(:,j));
                    temp = feval(obj.BaseFunction, distant, obj.Constant);
                    A(i,j) = temp;
                    A(j,i) = temp;
                    deriv =  obj.rbfphi_deriv(obj.Nodes(:,i), obj.Nodes(:,j));
                    for k = 1: dim
                            A(i+k*n,j) = obj.BaseFunctionDeriv(deriv(k),obj.Constant);
                            A(j+k*n,i) = A(i+k*n,j);
                    end
                end
                A(i,i) = A(i,i) - obj.Smooth;
            end
          
        end
        
        function [A] = AssembleGradBonus(obj)
           [dim, n] = size(obj.Nodes);
            A = zeros(n*(dim+1),n);
            for i = 1:n
                for j = 1:i
                    distant = norm(obj.Nodes(:,i) - obj.Nodes(:,j));
                    temp = feval(obj.BaseFunction, distant, obj.Constant);
                    A(i,j) = temp;
                    A(j,i) = temp;
                    deriv =  obj.rbfphi_deriv(obj.Nodes(:,i), obj.Nodes(:,j));
                    for k = 1: dim
                            A(i+k*n,j) = obj.BaseFunctionDeriv(deriv(k),obj.Constant);
                            A(j+k*n,i) = A(i+k*n,j);
                    end
                end
                A(i,i) = A(i,i) - obj.Smooth;
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
        function deriv =rbfphi_deriv(obj, xFixed, x)
            deriv = ones(size(x,1));
            norm = distant(x, xFixed);
            for j= 1:size(x,1)
                defiv(j) = (x(j) - xFixed(j))/norm;
            end
        end
        
        function u=rbfphi_linear_deriv(obj, r, const)
            u=1; end
        function u=rbfphi_cubic_deriv(obj, r, const)
            u = 2*r.*r;end
        function u=rbfphi_gaussian_deriv(obj, r, const)
            u = exp(-0.5*r.*r/(const*const));end
        function u=rbfphi_multiquadrics_deriv(obj, r, const)
            u=sqrt(1+r.*r/(const*const));end
        function u=rbfphi_thinplate_deriv(obj, r, const)
            u = r.*r.*log(r+1);end
        function u=rbfphi_cubicspline_deriv(obj, r, const)
            u= (r.*r + const*const).^(3/2);end
    end
end


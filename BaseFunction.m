classdef BaseFunction
    
    properties
        Func,
        Deriv
    end
    
    methods
        function obj = BaseFunction(name)
        switch lower(name)
                  case 'linear'          
                    obj.Func   = @obj.rbfphi_linear;
                    obj.Deriv = @obj.rbfphi_linear_deriv;
                  case 'cubic'
                    obj.Func   = @obj.rbfphi_cubic;
                    obj.Deriv =  @obj.rbfphi_cubic_deriv;
                  case 'multiquadric'
                    obj.Func   = @obj.rbfphi_multiquadrics;
                    obj.Deriv = @obj.rbfphi_multiquadrics_deriv;
                  case 'thinplate'
                    obj.Func   = @obj.rbfphi_thinplate;
                    obj.Deriv =  @obj.rbfphi_thinplate_deriv;
                  case 'gaussian'
                    obj.Func   = @obj.rbfphi_gaussian;
                    obj.Deriv = @obj.rbfphi_gaussian_deriv;
                  case 'cubicspline'
                    obj.Func   = @obj.rbfphi_cubicspline;
                    obj.Deriv = @obj.rbfphi_cubicspline_deriv;
                otherwise
                    obj.Func   = @obj.rbfphi_linear;
                    obj.Deriv =  @obj.rbfphi_linear_deriv;
            end    
        end
        function deriv = normDeriv(obj, xFixed, x)
            norm = distant(x, xFixed);
            for j= 1:size(x,1)
                deriv(j) = (xFixed(j) - x(j))/norm;
            end
        end
    end
    methods(Access = private)
        
        % Radial Base Functions
        %*************************************
         function u=rbfphi_linear(obj, r, const)
            u=r; end
        function u=rbfphi_cubic(obj, r, const)
            u=r.*r.*r;end
        function u=rbfphi_gaussian(obj, r, const)
            u=exp(-0.5*r.*r/(const*const));end
        function u=rbfphi_multiquadrics(obj, r, const)
            u=sqrt(1+r.*r/(const*const));end
        function u = rbfphi_thinplate(obj, r, const)
            u = r.*r.*log(r+1);end
        function u=rbfphi_cubicspline(obj, r, const)
            u= (r.*r + const*const).^(3/2);end
        
        % Radial Base Functions Derivatives
        
        
        function u=rbfphi_linear_deriv(obj, r, const)
            u=const; end
        function u=rbfphi_cubic_deriv(obj, r, const)
            u = 2*r.*r;end
        function u=rbfphi_gaussian_deriv(obj, r, const)
            u = -r/(const*const)*exp(-0.5*r.*r/(const*const));end
        function u=rbfphi_multiquadrics_deriv(obj, r, const)
            u= r/(sqrt(1+r.*r/(const*const))*const*const);end
        function u=rbfphi_thinplate_deriv(obj, r, const)
            u = 2*r*log(r+1) + (r*r)/(r+1);end
        function u=rbfphi_cubicspline_deriv(obj, r, const)
            u= 3*r* (r.*r + const*const).^(1/2);end
    end
    
end


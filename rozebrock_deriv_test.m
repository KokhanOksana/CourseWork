clear
count = 100; dim = 2;
a = -1.5; b = 1.5;
x0 = [2 2];
testFunc = Rozenbrock();

design = lhsdesign(count,dim);
x(1,:) = a + (b-a)* design(:,1); 
x(2,:) = a + (b-a)* design(:,2);
z = testFunc.Func(x);
derives = testFunc.Deriv(x);

ti = a:.05:b; 
[X,Y] = meshgrid(ti,ti);
F = testFunc.FuncMesh(X,Y);

subplot(1,3,1); mesh(X,Y,F); title('Real'); hold on;
%[x_minr, f_x_minr,ex_flagr, outputr] = fmincon(@testFunc.Func,[-1  2],[],[],[],[],[0 0],[1 1]);
%plot3(x_minr(1),x_minr(2),f_x_minr,'*b') ;hold on;

%RBF interpolation
    %varargin(1) - BaseFunctionName
    %varargin(2) - Constant
    %varargin(3) - Smooth
    %varargin(4) - DerivValues
rbf = RBF(x, z','cubic',0,0,derives);
Z = rbf.Interpolate([X(:)'; Y(:)']);
Z = reshape(Z, size(X));

subplot(1,3,2); mesh(X,Y,Z), hold
plot3(x(1,:), x(2,:), z,'*r'), hold off; title('RBF interpolation'); hold on; 
subplot(1,3,3); mesh(X,Y,abs(Z - F));  title('Error');

%figure;
%min_opt = optimset('TolX',1e-10);
%[x_minir, f_minir] = f_region_min(rozen,x0,25,[a a],[b b],min_opt)




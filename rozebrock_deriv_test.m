clear
figure;

count = 100; dim = 2;
aditionalCount = count/4;
a = -2; b = 2;
x0 = [2 2];
testFunc = Rozenbrock();

design = lhsdesign(count + aditionalCount,dim);
xAditional(1,:) = a + (b-a)* design(count+1:count+aditionalCount,1); 
xAditional(2,:) = a + (b-a)* design(count+1:count+aditionalCount,2);
zAditional = testFunc.Func(xAditional);
x(1,:) = a + (b-a)* design(1:count,1); 
x(2,:) = a + (b-a)* design(1:count,2);
z = testFunc.Func(x);
derives = testFunc.Deriv(x);

ti = a:.05:b; 
[X,Y] = meshgrid(ti,ti);
F = testFunc.FuncMesh(X,Y);

subplot(1,3,1); mesh(X,Y,F); title('Real'); hold on;
% [x_minr, f_x_minr,ex_flagr, outputr] = fmincon(@testFunc.Func,[-1  2],[],[],[],[],[0 0],[1 1]);
% plot3(x_minr(1),x_minr(2),f_x_minr,'*b') ;hold on;

%RBF interpolation
rbf = RBF(x, z', 'thinplate');
Z = rbf.Interpolate([X(:)'; Y(:)']);
Z = reshape(Z, size(X));

 subplot(1,3,2); mesh(X,Y,Z), hold
 plot3(x(1,:), x(2,:), z,'*r'), hold off; title('RBF interpolation'); hold on; 
 plot3(xAditional(1,:), xAditional(2,:), zAditional,'*b');hold on;
 subplot(1,3,3); mesh(X,Y,abs(Z - F));  title('Error');




clear
count = 100; n = 2;
a = -1.5; b = 1.5;
x0 = [2 2];
rozen = @(x)100*(x(2) - x(1).^2).^2 + (1-x(1)).^2; 
rozen2 = @(x,y)100*(y - x.^2).^2 + (1-x).^2;


X = lhsdesign(count,n);
x = a + (b-a)* X(:,1); 
y = a + (b-a)* X(:,2);
z = rozen2(x,y);

ti = a:.05:b; 
[XI,YI] = meshgrid(ti,ti);
FI = rozen2(XI, YI);

subplot(1,3,1); mesh(XI,YI,FI); title('Real'); hold on;
[x_minr, f_x_minr,ex_flagr, outputr] = fmincon(rozen,[-1  2],[],[],[],[],[0 0],[1 1]);
plot3(x_minr(1),x_minr(2),f_x_minr,'*b') ;hold on;

%RBF interpolation
%op=rbf_create([x'; y'], z','RBFFunction', 'thinplate');
rbf = RBF([x'; y'], z','thinplate');
ZI = rbf.Interpolate([XI(:)'; YI(:)']);
ZI = reshape(ZI, size(XI));

subplot(1,3,2); mesh(XI,YI,ZI), hold
plot3(x,y,z,'*r'), hold off; title('RBF interpolation'); hold on; 
subplot(1,3,3); mesh(XI,YI,abs(ZI - FI));  title('Error');
%
figure;
%min_opt = optimset('TolX',1e-10);
%[x_minir, f_minir] = f_region_min(rozen,x0,25,[a a],[b b],min_opt)




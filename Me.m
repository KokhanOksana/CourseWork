function [ x_min_interp, f_x_min_interp] = f_region_min( fun, x0, lb0, ub0 )
global dim d eps node_count k funccount  x_min_interp_prev ;
d =  max( (ub0 - lb0)/5);   dim = size(x0,2);   
node_count = 10;    eps = 0.0001;   k = 0;  funccount = node_count + 1;

x_min_interp = x0';   x_min_interp_prev = x_min_interp + 2*eps;
f_x_min_interp = fun(x_min_interp) + 2*eps;
lb = x_min_interp' - d*ones(1,dim);
ub = x_min_interp' + d*ones(1,dim);

%початкові точки експерименту    
[lb, ub] = bounds(x_min_interp, lb0, ub0) ;  
x = lg_design(lb, ub, node_count);

%початкові значення 
for j = 1:node_count
    z(j) = fun(x(j,:)');
end

dg = inf;
df = inf;

while( dg > eps)
    k = k + 1;

    %задання меж області
    lb_prev = lb;   ub_prev = ub;
    [lb, ub] = bounds(x_min_interp, lb, ub, lb0, ub0);

    %вибір наступних точок
    x_prev = x; z_prev = z;
    [x,z] = values(lb, ub, lb_prev, ub_prev, x_prev, z_prev, fun);
    
    %обрахування нового мінімуму
    x_min_interp_prev = x_min_interp;
    f_x_min_interp_prev = f_x_min_interp;
    [x_min_interp, f_x_min_interp] = f_sub_region_min( x_min_interp', x, z', lb,ub );
    
    %половина розміру сторонипрямокутника - області
    d = distant(x_min_interp, x_min_interp_prev);
    dg = abs((f_x_min_interp - f_x_min_interp_prev) / d); 
    df = abs(f_x_min_interp - f_x_min_interp_prev);
    if(d < max( (ub0 - lb0)/40))
        d = max( (ub0 - lb0)/40);
    end
    
    if(k > d*800)
        display('break because of step abuse');
        break;
    end;
    
    %зупинка , якщо наступна точка дуже близько до попередньої
    if(max(abs( x_min_interp_prev - x_min_interp))< 0.000001)
        display('break because of points closeness');
        break;
    end
end;
plot(x_min_interp(1),x_min_interp(2),'*b');hold on;
funccount
end

function [ x_min_interp, f_x_min_interp] = f_sub_region_min( x0, X, Z, lb, ub)
    %побудова RBF
    rbf = RBF(X', Z, 'thinplate'); 
    min_opt = optimset('Display','off');
    
    %знаходження мінімуму інтепполюючої функції
    [x_min_interp, f_x_min_interp] = fmincon(@rbf.Interpolate, x0',[],[],[],[], lb,ub,[], min_opt);
end

%обрізання області, якщо межі виходять за початкові
function [lb, ub, bound_size] = correct_bouns(lb, ub, lb0, ub0)
    global dim;    
    for l=1:dim
        if(lb(l) < lb0(l))
            lb(l) = lb0(l);
        end
        if(ub(l) > ub0(l))
            ub(l) = ub0(l);
        end
        b = size([lb(l):0.01:ub(l)]');
        bound_size(l) = b(1);
    end;
end

%обрахування меж області навколо точки
function [lb, ub] = bounds(x_min_interp, lb0, ub0)
global d dim;
    lb = x_min_interp' - d*ones(1,dim);
    ub = x_min_interp' + d*ones(1,dim);
    [lb, ub, bound_size ] = correct_bouns(lb, ub, lb0, ub0);
    
    % виведення на графіку меж
    plot(x_min_interp(1),x_min_interp(2),'*r');hold on;
    
    plot( [lb(1):0.01:ub(1)], ub(2)*ones(1,bound_size(1)) ,'-g');hold on;
    plot( [lb(1):0.01:ub(1)], lb(2)*ones(1,bound_size(1)) ,'-g');hold on;
    plot( lb(1)*ones(1,bound_size(2)), [lb(2):0.01:ub(2)] ,'-g');hold on;
    plot( ub(1)*ones(1,bound_size(2)), [lb(2):0.01:ub(2)] ,'-g');hold on;
end

%план ескперименту і переведення (0,1)->(a,b)
function x = lg_design(lb, ub, node_count)
global dim; 
    x_design = lhsdesign(node_count,dim);
    for j = 1:node_count
        for p = 1:dim
            x(j,p) = lb(p) + (ub(p)-lb(p))* x_design(j,p);
        end
    end
end

%відбір точок з попередньої облсті і побудова нових, обчаслення значення
%цільової функції в них
%алгоритм описано в курсовій
function [x, z] = values(lb, ub, lb_prev, ub_prev, x_prev, z_prev, fun)
global funccount node_count dim;
x_prob = lg_design(lb, ub, node_count);    
count = 0;
clear x;
for l = 1:size(x_prev,1)
    take_xprev = 1;
    take_xdesigm = 0;
            %plot(x_prev(l,1), x_prev(l,2),'.y')
    if(l <= size(x_prob,1)) 
            %plot(x_prob(l,1), x_prob(l,2),'.m') 
    end
    for k = 1:dim
        if((x_prev(l,k) < max(lb(k),lb_prev(k))) || (x_prev(l,k) > min(ub(k),ub_prev(k))))
            take_xprev = 0;
        end
        if(l <= size(x_prob,1)) 
            if((x_prob(l,k) > ub_prev(k)) || (x_prob(l,k) < lb_prev(k)))
                take_xdesigm = 1;
            end
        end
    end
    if(take_xprev == 1)
        count = count + 1;
        x(count,:) = x_prev(l,:);
        z(count) = z_prev(l);
            %plot(x(count,1), x(count,2), '.g'); hold on
    end
    if(take_xdesigm == 1)
        count = count + 1;
        x(count,:) = x_prob(l,:);
        z(count) = fun(x_prob(l,:)');
        funccount = funccount+1;
            %plot(x(count,1), x(count,2), '.b'); hold on
    end
end;

if(norm(x_min_interp - x_min_interp_prev) < 0.05)
    x_prob = lg_design(lb, ub, 4);
    for l = 1:size(x_prob)
        count = count+1;
        x(count,:) = x_prob(l,:);
        z(count) = fun(x_prob(l,:)');
        funccount = funccount+1;
            %plot(x(count,1), x(count,2), '.b'); hold on
    end
end
%добирання точок , якщо не вистарчає до плану експерименту
if(count < node_count)
    x_prob = lg_design(lb, ub, node_count - count);
    for l = 1:size(x_prob)
        count = count+1;
        x(count,:) = x_prob(l,:);
        z(count) = fun(x_prob(l,:)');
        funccount = funccount+1;
            %plot(x(count,1), x(count,2), '.b'); hold on
    end
end
end

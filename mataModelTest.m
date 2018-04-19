a = -2; b = 2;
x0 = [2 2];
testFunc = Rozenbrock();
[x_minir, f_minir] = f_region_min(@testFunc.Func,x0,[a a],[b b])


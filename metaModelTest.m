eps = 0.000001;
nodeCount = 10;
x0 = [1.5 0.5]';
global dim callCount;
dim = size(x0,1);
callCount = 0; 
bounds = Bounds([0 0],[2 2]);
testFunc = Rozenbrock();

metaModel = MetaModelDeriv(testFunc, eps, nodeCount);

metaModel.Min(x0,bounds).ShowResult();

%% TEST CASE: DISTORTION OPERATOR AS A 1D MARGINALIZATION
T = sptensor(tensor(@ones, [10 10 10]));
v = ssit.FspVector(T);

f1 = @(x1, y1) (y1 == x1);
f2 = @(x2, y2) 1;
f3 = @(x3, y3) 1;

cpds = {f1, f2, f3};
yranges = {0:9, [1], [1]};

perms = {[1,2,3], [3,2,1], [2,3,1], [1,3,2], [3,1,2], [2,1,3]};
for i = 1:length(perms)
    perm = perms{i};
    C = ssit.pdo.TensorProductDistortionOperator(cpds(perm), yranges(perm));
    w = C.computeObservationDist(v);
    assert(all(double(squeeze(w.data)) == 10*10*ones(10, 1)));
end

%% TEST CASE: DISTORTION OPERATOR AS A 2D MARGINALIZATION
T = sptensor(tensor(@ones, [100 100 100]));
v = ssit.FspVector(T);

f1 = @(y1, x1) (y1 == x1);
f2 = @(y2, x2) (y2 == x2);
f3 = @(y3, x3) 1;

cpds = {f1, f2, f3};
yranges = {0:99, 0:99, [1]};

perms = {[1,2,3], [3,2,1], [2,3,1], [1,3,2], [3,1,2], [2,1,3]};
for i = 1:length(perms)
    perm = perms{i};
    C = ssit.pdo.TensorProductDistortionOperator(cpds(perm), yranges(perm));
    w = C.computeObservationDist(v);
    assert(all(all(double(squeeze(w.data)) == 100*ones(100, 100))));
end

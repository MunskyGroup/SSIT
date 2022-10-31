%% TEST CASE: ADDING SPARSE VECTORS

X1 = [
    0 0 1;
    1 2 3;
    0 2 1;
]';

X2 = [
    0 0 1;
    0 0 2;
    1 2 3;
    0 0 4;
    ]';

p1 = [1 1 1]';

p2 = [1 0 0 3]';

X_add = [
    0 0 1;
    1 2 3;
    0 2 1;
    0 0 2;
    0 0 4;
    ]';
p_add = [2 1 1 0 3]';

[X_new, p_new] = ssit.math.distrAdd(X1, p1, X2, p2);

assert(isequal(X_new, X_add));
assert(isequal(p_new, p_add));
function m = smolyakLevelSize(level)
% m = ssit.pgf.smolyakLevelSize(level) number of 1-D Chebyshev-Gauss-Lobatto
% nodes at a (1-based) Smolyak level: 1 for level 1, 2^(level-1)+1 otherwise.
m = 2.^(level-1) + 1;
m(level == 1) = 1;
end

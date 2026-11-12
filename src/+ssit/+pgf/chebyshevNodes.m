function x = chebyshevNodes(m)
% x = ssit.pgf.chebyshevNodes(m) returns the m Chebyshev-Gauss-Lobatto
% points (Chebyshev extrema) on [-1,1] as a column vector, in ascending
% order: x(i) = -cos(pi*(i-1)/(m-1)).  m = 1 returns 0.
%
% The construction mirrors ChebyshevApprox.chebyshev_extrema (used by
% SmolyakApprox.jl) so that the nested levels of the sparse grid contain
% bit-identical nodes: the argument (i-1)/(m-1) of node i at level m equals
% the argument of node 2i-1 at level 2m-1 exactly, and the middle node is 0.
arguments
    m (1,1) double {mustBeInteger, mustBePositive}
end
x = zeros(m,1);
for i = 1:floor(m/2)
    xi = -cospi((i-1)/(m-1));
    x(i) = x(i) + xi;
    x(m-i+1) = x(m-i+1) - xi;
end
end

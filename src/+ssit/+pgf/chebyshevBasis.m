function P = chebyshevBasis(order, x, derivOrder)
% P = ssit.pgf.chebyshevBasis(order, x, derivOrder) returns the
% numel(x) x (order+1) matrix whose column n+1 holds the derivOrder-th
% derivative of the Chebyshev polynomial T_n evaluated at x (in [-1,1]).
%   derivOrder = 0 (default), 1 or 2.
% Uses the three-term recurrence T_{n+1} = 2 x T_n - T_{n-1} and its first
% and second derivatives (as in ChebyshevApprox.jl / Edward's utils.jl).
arguments
    order (1,1) double {mustBeInteger, mustBeNonnegative}
    x (:,1) double
    derivOrder (1,1) double {mustBeMember(derivOrder,[0 1 2])} = 0
end
N = numel(x);
T = zeros(N, order+1);
dT = zeros(N, order+1);
ddT = zeros(N, order+1);
T(:,1) = 1;
if order >= 1
    T(:,2) = x;
    dT(:,2) = 1;
end
for n = 2:order
    T(:,n+1) = 2*x.*T(:,n) - T(:,n-1);
    dT(:,n+1) = 2*T(:,n) + 2*x.*dT(:,n) - dT(:,n-1);
    ddT(:,n+1) = 4*dT(:,n) + 2*x.*ddT(:,n) - ddT(:,n-1);
end
switch derivOrder
    case 0
        P = T;
    case 1
        P = dT;
    case 2
        P = ddT;
end
end

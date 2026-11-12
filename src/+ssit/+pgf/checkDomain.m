function dom = checkDomain(domain, d)
% dom = ssit.pgf.checkDomain(domain, d) validates a sparse-grid domain and
% returns it as a 2 x d matrix with row 1 = upper bounds and row 2 = lower
% bounds (the SmolyakApprox.jl convention).  A column given as
% [lower; upper] is swapped, so either orientation is accepted.
arguments
    domain (2,:) double
    d (1,1) double {mustBeInteger, mustBePositive}
end
if size(domain,2) ~= d
    error('ssit:pgf:domainSize', ...
        'The domain has %d columns but the grid has %d dimensions.', size(domain,2), d);
end
dom = domain;
swap = dom(1,:) < dom(2,:);
dom(:,swap) = dom([2 1],swap);
end

function L = LogLikelyhood(P_func, Q_func, x)
    % LogLikelihood computes the log likelihood between two probability functions P and Q
    % Input:
    %   P_func - a function handle representing the probability function P
    %   Q_func - a function handle representing the probability function Q
    %   x - a vector of points at which to evaluate the probability functions
    % Output:
    %   L - the log likelihood value

    % Evaluate the probability functions at the given points
    P = P_func(x);
    Q = Q_func(x);

    P = P / sum(P);
    Q = Q / sum(Q);

    P = P + eps;
    Q = Q + eps;
    
    % Check if the evaluated distributions are valid
    if any(P < 0) || any(Q < 0)
        error('Probability functions must return non-negative values.');
    end
    % Calculate the log likelihood
    L = sum(P .* log(Q));
end
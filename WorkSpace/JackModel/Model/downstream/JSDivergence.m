function d = JSDivergence(P_func, Q_func, x_range)
    P = P_func(x_range);
    Q = Q_func(x_range);

    P = P / sum(P);
    Q = Q / sum(Q);


    P = P + eps;
    Q = Q + eps;

    M = 0.5 * (P + Q);
    d = 0.5 * sum(P .* log2(P ./ M)) + 0.5 * sum(Q .* log2(Q ./ M));
end
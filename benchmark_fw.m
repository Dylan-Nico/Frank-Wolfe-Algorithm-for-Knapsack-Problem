n = 10
eps = 1e-6

for k = 0:9
    seed = k;              % seed diverso a ogni run
    rng(seed,'twister');   % riproducibilità

    [Q, q, l, u, a, b, x0] = generate_variables(n, seed);
    Q = (Q + Q')/2;
    frank_wolfe(Q, q, x0, a, b, l, u, eps)
    frank_wolfe_away(Q, q, x0, a, b, l, u, eps)
end

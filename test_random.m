n = randi([2, 100]);
[Q, q, l, u, a, b, x0] = generate_variables(n);
eps = 1e-8;

[x_FW, f_FW] = frank_wolfe(Q, q, x0, a, b, l, u, eps, "normal");
[x_FW, f_FW] = frank_wolfe(Q, q, x0, a, b, l, u, eps, "away");
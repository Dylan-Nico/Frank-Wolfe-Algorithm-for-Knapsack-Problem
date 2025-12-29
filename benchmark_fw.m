n = 10
eps = 1e-6

for k = 0:5
    seed = k;              % seed diverso a ogni run
    rng(seed,'twister');   % riproducibilità

    for type = ["Interior", "Boundary"]

        [Q, q, l, u, a, b, x0, x_star] = generate_variables(n, seed, type);
        [x_FW, f_FW, f_star, gaps_FW, pe_FW] = frank_wolfe(Q, q, x0, a, b, l, u, eps, x_star);
        [x_AW, f_AW, f_star2, gaps_AW, pe_AW] = frank_wolfe_away(Q, q, x0, a, b, l, u, eps, x_star);
        comparison_plot(gaps_FW, gaps_AW, pe_FW, pe_AW, type, n);
        
    end
end

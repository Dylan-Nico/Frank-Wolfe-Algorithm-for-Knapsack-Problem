function [Q, q, l, u, a, b, x0, x_star] = generate_variables(n, seed, type)

    rng(seed)
    
    % -------------------------
    % Domain
    % -------------------------
    l = zeros(n,1);
    u = ones(n,1);
    a = rand(n,1) + 0.1;
    
    
    % -------------------------
    % Build SPD Q with controlled spectrum
    % -------------------------
    kappa = 1e2;
    lam_min = 1;
    lam_max = kappa * lam_min;
    
    lam = logspace(log10(lam_min), log10(lam_max), n)';
    [U,~] = qr(randn(n),0);
    Q = U * diag(lam) * U';
    Q = 0.5 * (Q + Q');
    if isequal(type, "Interior")
        x_star = 0.3 + 0.4*rand(n,1);
    else 
        if isequal(type, "Boundary")
            x_star = ones(n,1);
            x_star(1) = 0.8;
        end
    end
    q = -2*Q*x_star;
    b = 0.1;
    x0 = rand(n,1);
end

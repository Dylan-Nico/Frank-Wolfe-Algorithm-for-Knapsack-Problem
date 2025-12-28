function generate_test_case(n, type)
% Create and RUN a test instance for Frank-Wolfe.
%
% Usage:
%   generate_test_case(n, "opt_on_vertex")
%   generate_test_case(n, "interior")
%   generate_test_case(n, "box_boundary")
%   generate_test_case(n, "active_linear")
%   generate_test_case(n, "interior_zero")
%   generate_test_case(n, "conditioning_vs_iterations")

% set seed for reproducibility 
rng(43);

fprintf("--------------------------------------------------\n");
fprintf(" Running test (%s), n = %d\n", type, n);
fprintf("--------------------------------------------------\n");



% -------------------------
% Domain
% -------------------------
l = zeros(n,1);
u = ones(n,1);
a = rand(n,1) + 0.1;

% -------------------------
% conditioning vs iterations (on interior case only)
% -------------------------
if strcmpi(type, "conditioning_vs_iterations")

    % conditioning number
    kappa_list = [1e2 1e3 1e4 1e5 1e6];
    
    % eigenvalues
    lam_min = 1;
    eps = 1e-3;

    % Fix x_star and x0 across kappas to make comparison fair
    x_star = 0.3 + 0.4*rand(n,1);
    x0     = rand(n,1);
    b      = 0.1;  % linear constraint inactive

    iters = zeros(length(kappa_list),1);
    final_gaps = zeros(length(kappa_list),1);

    for i = 1:length(kappa_list)
        kappa = kappa_list(i);
        lam_max = kappa * lam_min;

        lam = logspace(log10(lam_min), log10(lam_max), n)';
        [U,~] = qr(randn(n),0);
        Q = U * diag(lam) * U';
        Q = 0.5 * (Q + Q');

        q = -2*Q*x_star;

        [x_FW, f_FW, f_star, gaps, primal_errors] = frank_wolfe(Q, q, x0, a, b, l, u, eps);

        iters(i) = length(gaps);

        % save gaps for eaxh run
        final_gaps(i) = gaps(end);

        fprintf("k=%.2e | iters=%d | f_FW=%.6f | f_star=%.6f\n", kappa, iters(i), f_FW, f_star);
    end

    % Plot
    figure; loglog(kappa_list, iters, '-o', 'LineWidth', 1.5);
    grid on;
    xlabel('k'); ylabel('iterations');
    title('Iterations vs conditioning (interior)');

    % plot reached gaps
    figure;
    loglog(kappa_list, final_gaps, '-o', 'LineWidth', 1.5);
    grid on;

    yline(eps, '--', 'eps', 'LineWidth', 1.2);

    xlabel('k');
    ylabel('final duality gap');
    title('Final duality gap vs conditioning (interior)');


    fprintf("--------------------------------------------------\n");
    return;
end

% -------------------------
% Choose conditioning 
% -------------------------
switch lower(type)
    case {'interior', 'box_boundary', 'active_linear','opt_on_vertex', 'interior_zero'}
        kappa = 1e2;
        cond_label = "well-conditioned";

    otherwise
        error("Test type not recognized.");
end

% -------------------------
% Build SPD Q with controlled spectrum
% -------------------------
lam_min = 1;
lam_max = kappa * lam_min;

lam = logspace(log10(lam_min), log10(lam_max), n)';
[U,~] = qr(randn(n),0);
Q = U * diag(lam) * U';
Q = 0.5 * (Q + Q');

% -------------------------
% Select geometry
% -------------------------
switch lower(type)

    case {'interior'}
        description = sprintf("Interior optimum, (%s Q).", cond_label);
        x_star = 0.3 + 0.4*rand(n,1);
        q = -2*Q*x_star;
        b = 0.1;

    case {'box_boundary'}
        description = sprintf("Optimum on box boundary (%s Q).", cond_label);
        x_star = ones(n,1);
        x_star(1) = 0.8;
        q = -2*Q*x_star;
        b = 0.1;

    case {'active_linear'}
        description = sprintf("Optimum on active linear constraint (%s Q).", cond_label);
        x_star = rand(n,1);
        b = a' * x_star;
        q = -2*Q*x_star;
        
    case {'opt_on_vertex'}
            description = sprintf("Optimum on vertex (%s Q).", cond_label);
            x_star = u; 
            q = -2*Q*x_star;
            b = 0.1;

    case {'interior_zero'}
        description = sprintf("Interior optimum with PSD Q (zero eigenvalues).");
        x_star = 0.3 + 0.4*rand(n,1);

        % Some eigenvalues equal to 1e-8
        lam(1:floor(n/4)) = 1e-8;               
        Q = U * diag(lam) * U';
        Q = 0.5 * (Q + Q');

        q = -2*Q*x_star;
        b = 0.1;

  
end

% -------------------------
% Starting point
% -------------------------
x0 = rand(n,1);

% -------------------------
% Precision
% -------------------------
eps = 1e-6;


% Run Frank–Wolfe (standard)
% -------------------------
[x_FW, f_FW, f_star, gaps_FW, pe_FW] = frank_wolfe(Q, q, x0, a, b, l, u, eps, x_star);

% -------------------------
% Run Frank–Wolfe (away)
% -------------------------
[x_AW, f_AW, f_star2, gaps_AW, pe_AW] = frank_wolfe_away(Q, q, x0, a, b, l, u, eps, x_star);


% Comparison plot
comparison_plot(gaps_FW, gaps_AW, pe_FW, pe_AW, type, n);


% -------------------------
% Print
% -------------------------
e = eig(Q);
fprintf("Q conditioning: k ≈ %.2e (λ_min=%.2e, λ_max=%.2e)\n", ...
        max(e)/min(e), min(e), max(e));

fprintf("Description: %s\n", description);
fprintf("f(x_FW)     = %.6f\n", f_FW);
fprintf("f(x_AW)     = %.6f\n", f_AW);
fprintf("f_star      = %.6f\n", f_star);
fprintf("--------------------------------------------------\n");

end



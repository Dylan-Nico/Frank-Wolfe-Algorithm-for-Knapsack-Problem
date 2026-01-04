function [x, fval, f_star, gaps, primal_errors] = frank_wolfe_away(Q,q,x0,a,b,l,u,eps,x_star_true)
%% Invariant:
% 1) keep the starting point x_0 feasible
% 2) the starting point must be a convex comb of vertex: x_0 = V*lambda, lambda >= 0 , sum_i(lambda_i)=1
%
% Away-Step Frank-Wolfe algorithm 
%
% Input: Q n*n Positive Semidefinite matrix
%        q n vector
%        x0 starting point (could be feasible or not, a procedure will
%        force it) --> thi version starts with closest point x0 feasible
%        for away step.
%        a,b for linear constraint
%        l,u margins of the box
%        eps precision - 10^-6 default
%
% Output: x solution
%         fval value of function in the solution point
%         f_star true value of function given by the Oracle

% Initialize variables

x = x0;
n = length(x0);
max_iter = 100000;

iterates = x;
gap = Inf;
gaps = [];
primal_errors = [];
Results = table([], [], [], [], [], [], 'VariableNames', {'Iter','gap','alpha','a * x','StepNorm', 'PrimalError'});


% check if starting point is feasible, if not, force it 
if ~check_feasible(x, a, b, l, u)
    fprintf("Error: starting point no feasible, force it to be feasible\n");
    
    % put in the center of the box
    x = (l + u) / 2;

    % if second constraint not verified
    if a' * x < b
        % Compute difference and move
        t = (b - a' * x) / (a' * a);
        x = x + t * a;

        % re-put in the box if necessary
        x = min(max(x, l), u);
    end
else
    fprintf("Starting point is feasible\n");
end

% initialize AFW with 2 feasible vertixes (no full decomposition)

% gradient in x0
g0 = 2*Q*x + q;

% pick two feasible vertixes
v1 = solveLP(g0,  a, b, l, u);   
v2 = solveLP(-g0, a, b, l, u);   

% build  2-vertex active set
V = [v1, v2];

% if v1 == v2, start at that single vertex
diffv = v1 - v2;
den = diffv' * diffv;

if den <= 1e-30
    lambda = 1;
    V = v1;
    x = v1;
else
    % project x in the segment [v2, v1] to get weights
    alpha = ((x - v2)' * diffv) / den;
    alpha = min(1, max(0, alpha));   

    lambda = [alpha; 1 - alpha];
    x = V * lambda;                  % decomposition
end


% Compute f_star: use true x_star if provided, otherwise call Oracle
if nargin >= 9 && ~isempty(x_star_true)
    f_star = x_star_true' * Q * x_star_true + q' * x_star_true;
else
    [~, f_star] = Oracle(Q, q, a, b, l, u);
end

k = 1;
tol_drop = 1e-12;

while (k < max_iter && gap > eps)

    % gradient
    g = 2*Q*x + q;

    % FW vertex 
    s = solveLP(g, a, b, l, u);
    d_fw = s - x;
    gap_fw = -g' * d_fw;   

    % Away vertex 
    % choose vertex v in V maximizing <g, v>  
    inner = V' * g;                
    [~, idxA] = max(inner);
    v = V(:, idxA);
    d_away = x - v;
    gap_away = -g' * d_away;      

    % Choose d
    if gap_fw >= gap_away  
        d = d_fw;
        alpha_max = 1;
        step_type = "FW";
    else 
        d = d_away;
        alpha_max = lambda(idxA) / (1 - lambda(idxA));
        step_type = "AWAY";
    end

    % use actual iteration gap
    gap_best = max(gap_fw,gap_away);

    % use still FW gap for certificare (bounded by theory)
    gap = gap_best;
    gaps(end+1) = gap;

    % ---- Exact line search: implicity usage of tomography (for quadratic function it's easy) ----
    % α = argmin f(x + α d) = phi(α) 
    num = g' * d;
    denom = 2*d' * Q * d;

    if denom > 0
        alpha_ls = -num / denom;
        alpha = min(alpha_max, max(0, alpha_ls));
    else
        alpha = alpha_max;
    end

    % step: update x
    x = x + alpha*d;

    % Keep the invariant true
    % Update active set weights 
    if step_type == "FW"
        % scale all existing weights by (1-alpha)
        lambda = (1 - alpha) * lambda;

        % check if vertex s already in V, if not, add new vertex s to
        % active set V
        found = false;
        for j = 1:size(V,2)
            if norm(V(:,j) - s) <= 1e-14 * max(1, norm(s))
                lambda(j) = lambda(j) + alpha;
                found = true;
                break;
            end
        end
        if ~found
            V = [V, s];
            lambda = [lambda; alpha];
        end

    else
        lambda = (1 + alpha) * lambda;
        lambda(idxA) = lambda(idxA) - alpha;

        % drop step for vertex in index idxA
        if lambda(idxA) <= tol_drop
            V(:, idxA) = [];
            lambda(idxA) = [];
        end
    end

    % numerical problem: clip tiny negatives and renormalize
    lambda(lambda < 0 & lambda > -1e-14) = 0;
    if isempty(lambda)
        V = solveLP(zeros(n,1), a, b, l, u);
        lambda = 1;
        x = V;
    else
        lambda = lambda / sum(lambda);
    end

    % evaluate f at current iteration and compute primal error
    f_x = x' * Q * x + q' * x;
    primal_error = abs(f_x - f_star) / max(1, abs(f_star));
    primal_errors(end+1) = primal_error;

    % save results
    iterNum   = k;
    gap_k     = gap;
    alpha_k   = alpha;
    aTx_k     = a' * x;
    stepnorm  = norm(alpha*d);
    pe_k      = primal_error;

    newRow = {iterNum, gap_k, alpha_k, aTx_k, stepnorm, pe_k};
    Results = [Results; newRow];

    iterates(:, end+1) = x;
    k = k + 1;
end

% save last iteration
gaps(end+1) = gap;
primal_errors(end+1) = primal_error;

fval = x' * Q * x + q' * x;

if n == 2
    FW_plot2D(Q, q, a, b, l, u, iterates);
end

format short e
disp(Results);
format short

if(check_feasible(x, a, b, l, u))
    fprintf("Solution point is feasible\n");
else
    fprintf("Solution point is not feasible\n");
end

plot_gaps(gaps, primal_errors);

end


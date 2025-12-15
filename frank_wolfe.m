function [x, fval, f_star, gaps, primal_errors] = frank_wolfe(Q,q,x0,a,b,l,u,eps) % invariant Sum_i(a_i*x_i)>=b

% Frank Wolfe algorithm
%
% Input: Q n*n Positive Semidefinite matrix
%        q n vector
%        x0 starting point (could be feasible or not, a procedure will force it)
%        a,b for linear constraint
%        l,u margins of the box
%        eps precision - 10^-3 default
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

% Compute true minimum with oracle (for primal error)
[x_star, f_star] = Oracle(Q, q, a, b, l, u);

k = 1;
while(k<max_iter && gap>eps)

    % gradient
    g = 2*Q*x + q;

    % Solve linear problem
    s = solveLP(g,a,b,l,u);
    
    % set direction
    d = s - x;
    
   % Compute duality gap and save it for the plot
    gap = -g' * d;
    gaps(end+1) = gap;


    % ---- Exact line search: implicity usage of tomography (for quadratic function it's easy) ----
    % α = argmin f(x + α d) = phi(α) 
    num = g' * d;
    denom = 2*d' * Q * d;

    if denom > 0
        alpha = min(1, max(0, -num / denom));
    else
        alpha = 1;
    end

    % step: update x
    x = x + alpha*d;

    % evaluate f at current iteration and compute primal error
    f_x = x' * Q * x + q' * x;
    primal_error = abs(f_x - f_star) / max(1, abs(f_x));
    primal_errors(end+1) = primal_error;

    % Save values
    iterNum   = k;
    gap_k     = gap;
    alpha_k   = alpha;
    aTx_k     = a' * x;
    stepnorm  = norm(alpha*d);
    pe_k      = primal_error;

    % Add values to table
    newRow = {iterNum, gap_k, alpha_k, aTx_k, stepnorm, pe_k};
    Results = [Results; newRow];
    
    % update next iterate
    iterates(:, end+1) = x;
    k = k+1;

end

% save gaps of last iteration
gaps(end+1) = gap;
primal_errors(end+1) = primal_error;

% evaluate f in last iteration
fval = x' * Q * x + q' * x;

% ---- Plot iterates and level sets is n=2 ----
if n == 2
    FW_plot2D(Q, q, a, b, l, u, iterates);
end

% print table
format short e
disp(Results);
format short

% check x ammissible
if(check_feasible(x, a, b, l, u))
    fprintf("Solution point is feasible\n");
else
    fprintf("Solution point is not feasible\n");
end


% plot gaps
plot_gaps(gaps, primal_errors);

end
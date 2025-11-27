function [x, fval] = frank_wolfe(Q,q,x0,a,b,l,u,eps) % invariant Sum_i(a_i*x_i)>=b

plot_tomo = false;
x = x0;
n = length(x0);
iterates = x;
gaps = [];

Results = table([], [], [], [], [], 'VariableNames', {'Iter','gap','alpha','a * x','StepNorm'});

% check if starting point is feasible
if ~check_feasible(x, a, b, l, u)
    fprintf("Error: starting point no feasible, choose one that respect the constraints\n");
    return;
end


for k = 0:5000
    g = 2*Q*x + q;

    s = solveLP(g,a,b,l,u);

    d = s - x;
    
    gap = -g' * d;
    % save gap for the plot
    gaps(end+1) = gap;

    if(gap <= eps)
        fprintf('Converged (gap <= eps)\n');
        break;
    end

    % ---- Line search esatta per quadratiche usando la tomografia implicitamente ----
    % α = argmin f(x + α d)
    num = g' * d;
    denom = 2*d' * Q * d;

    if denom > 0
        alpha = min(1, max(0, -num / denom));
    else
        alpha = 1;
    end

    if plot_tomo && mod(k,10)==0
        plot_tomography(Q, q, x, d, alpha);
    end

    x = x + alpha*d;

    % Valori salvati
    iterNum   = k;
    gap_k     = gap;
    alpha_k   = alpha;
    aTx_k     = a' * x;
    stepnorm  = norm(alpha*d);

    % Aggiungi alla tabella
    newRow = {iterNum, gap_k, alpha_k, aTx_k, stepnorm};
    Results = [Results; newRow];

    % === WARNING VICINO AL BOUND ===
    tolBound = 1e-8;
    warning_bound(tolBound,x,l,u) ;

    iterates(:, end+1) = x;

end

fval = x' * Q * x + q' * x;

% ---- Plot 2D delle iterazioni se n = 2 ----
if n == 2
    FW_plot2D(Q, q, a, b, l, u, iterates);
end

% stampo la tabella
disp(Results);

% check x ammissible
check_feasible(x, a, b, l, u);

% plot gap
plot_gap(gaps);

end
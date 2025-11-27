function [x, fval] = frank_wolfe(Q,q,x0,a,b,l,u,eps) % invariant Sum_i(a_i*x_i)>=b

x = x0;
n = length(x0);
iterates = x;

for k = 0:5000
    g = 2*Q*x + q;

    s = solveLP(g,a,b,l,u);

    d = s - x;
    
    gap = -g' * d;
    if(gap <= eps)
        break;
    end

    % ---- Line search esatta per quadratiche ----
    % α = argmin f(x + α d)
    num = g' * d;
    denom = 2*d' * Q * d;

    if denom > 0
        alpha = min(1, max(0, -num / denom));
    else
        alpha = 1;
    end

    x = x + alpha*d;
    iterates(:, end+1) = x;

end

fval = x' * Q * x + q' * x;

% ---- Plot 2D delle iterazioni se n = 2 ----
if n == 2
    FW_plot2D(Q, q, a, b, l, u, iterates);
end

end
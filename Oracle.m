function [x_star, f_star] = Oracle(Q, q, a, b, l, u)
% Oracle: find the optimal solution of the quadratic knapsack problem
% use quadprog
%
% min  x' Q x + q' x          
% s.t. a' x >= b
%      l <= x <= u

    n = length(q);

    % quadprog use 1/2 x' H x + f' x
    H = 2*Q;       % Because your f(x)=x'Qx + q'x, quadprog wants 1/2 x'Hx + f'x
    f = q;

    % Convert a'x >= b   →   -a'x <= -b
    Aineq = -a';
    bineq = -b;

    % Bounds
    lb = l;
    ub = u;

    opts = optimoptions('quadprog', 'Display', 'none');

    [x_star, f_star_qp] = quadprog(H, f, Aineq, bineq, [], [], lb, ub, [], opts);

    % Convert quadprog objective to final value f^*
    f_star = x_star' * Q * x_star + q' * x_star;
end

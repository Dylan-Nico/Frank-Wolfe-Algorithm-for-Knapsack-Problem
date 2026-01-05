function ok = check_feasible(x, a, b, l, u)
% Check if x is within the bounds defined by l and u and if linear
% constraint is valid
    tol = 1e-9;
    ok = (a' * x >= b - tol) && all(x >= l - tol) && all(x <= u + tol);
end

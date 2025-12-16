function ok = check_feasible(x, a, b, l, u)
    tol = 1e-9;
    ok = (a' * x >= b - tol) && all(x >= l - tol) && all(x <= u + tol);
    
    if ok
        disp("x feasible: true");
    else
        disp("x feasible: false");
    end
end

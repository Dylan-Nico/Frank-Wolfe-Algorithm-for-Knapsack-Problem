function warning_bound(tolBound, x, l, u)

if any(abs(x - l) < tolBound)
        fprintf('WARNING: touching LOWER bound\n');
end

if any(abs(x - u) < tolBound)
        fprintf('WARNING: touching UPPER bound\n');
end
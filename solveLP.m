function s = solveLP(g, a, b, l, u)
    n = length(g);

    % First: unconstrained box minimizer
    y = zeros(n,1);
    for i = 1:n
        if g(i) > 0
            y(i) = l(i);
        else
            y(i) = u(i);
        end
    end
    
    % If feasible, we are done
    if a'*y >= b
        s = y;
        return;
    end

    % Otherwise we must increase value a^T y until reaching b
    % We do so by moving along direction that increases a_i / g_i efficiency.
    
    % Compute "efficiency": how much a grows per increase of objective
    ratio = a ./ abs(g + 1e-12);

    % Sort variables by ratio (descending)
    [~, idx] = sort(ratio, 'descend');

    s = y;

    deficit = b - a'*s;

    for k = 1:n
        i = idx(k);
        if a(i) <= 0, continue; end

        maxInc = u(i) - s(i);
        if maxInc <= 0, continue; end

        gain = a(i)*maxInc;

        if gain >= deficit
            % only partial fill needed
            s(i) = s(i) + deficit / a(i);
            return;
        else
            % take all possible increase
            s(i) = u(i);
            deficit = deficit - gain;
        end
    end

    % should never reach here: if box is insufficient, problem is infeasible
end

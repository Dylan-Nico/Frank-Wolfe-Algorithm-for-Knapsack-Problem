function s = solveLP(g, a, b, l, u)
% solve knapsack problem
% Input: gradient g
%        a, b for linear constraint
%        l, u extremes of the box
% Output: s_k = argmin g^T*s, with s in the Domain

    n = length(g);

    % Unconstrained box minimizer
    % easy: if my line is increasing the minumum is l, if it's decreasing
    % the min is u
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
    % We move along direction that increases a_i / g_i in an efficient way
    
    % Compute "efficient": how much a grows per increase of objective
    ratio = a ./ abs(g + 1e-12);

    % Sort variables by ratio 
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

    % never reach here: if box is insufficient, problem is infeasible
end


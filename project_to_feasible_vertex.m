function [V, lambda] = project_to_feasible_vertex(x, a, b, l, u)
% Decompose a feasible point x into convex combination of vertices of:
%   { x : a' x >= b,  l <= x <= u }.
% Vertices have all variables at bounds except possibly one.
%
% Output:
%   V      n-by-m matrix of vertices (columns)
%   lambda m-by-1 convex weights, sum(lambda)=1

if nargin < 6 || isempty(tol), tol = 1e-12; end

tol = 1e-9;
n = length(x);
d = u - l;

free = (abs(d) > 0);
x_fixed = x; 

% Normalize y in [0,1] on free components
y = zeros(n,1);
y(free) = (x(free) - l(free)) ./ d(free);

% Clip numerics
y(y < 0 & y > -1e-12) = 0;
y(y > 1 & y < 1+1e-12) = 1;

% Check box feasibility 
if any(y(free) < -1e-8) || any(y(free) > 1+1e-8)
    error("x is not within [l,u] (up to tolerance).");
end

% Knapsack in y: c*y >= b
c = zeros(n,1);
c(free) = a(free) .* d(free);
B = b - a' * l;

% Check feasibility 
if c' * y + 1e-10 < B
    error("x is not feasible for a' x >= b (up to tolerance).");
end

% We will maintain a list of points y_k and weights w_k
Y = y;
W = 1;

while true
    % Find an index k of a point with >=2 fractional components
    found = false;
    k_split = -1;
    for k = 1:size(Y,2)
        yy = Y(:,k);
        frac = find( yy > tol & yy < 1 - tol & free ); % for free vars
        if numel(frac) >= 2
            found = true;
            k_split = k;
            break
        end
    end
    if ~found
        break
    end

    yy = Y(:,k_split);
    frac = find( yy > tol & yy < 1 - tol & free );
    i = frac(1);
    j = frac(2);

    % Need c_i, c_j > 0 for the "keep c'y constant" direction.
    % If some a_i or d_i is zero, c_i=0 -> it doesn't matter for knapsack,
    % just push it to bounds without affecting feasibility.
    if c(i) <= tol || c(j) <= tol
        t1 = yy(i);        
        t2 = 1 - yy(i);    
        y0 = yy; y0(i)=0;
        y1 = yy; y1(i)=1;
        theta = t2/(t1+t2); 
        w = W(k_split);

        % replace k_split with y0 and y1
        Y(:,k_split) = y0;
        W(k_split) = w*theta;
        Y = [Y, y1];
        W = [W; w*(1-theta)];
        continue
    end

    % Direction d keeps c' y constant:
    % y_i += t / c_i,  y_j -= t / c_j
    t_plus  = min( (1 - yy(i)) * c(i), yy(j) * c(j) );          % max t>0
    t_minus = min( yy(i) * c(i), (1 - yy(j)) * c(j) );          % max t<0 in magnitude

    if t_plus <= tol || t_minus <= tol
        % numerical degeneracy: snap one variable
        if yy(i) < 0.5
            yy(i)=0;
        else
            yy(i)=1;
        end
        Y(:,k_split)=yy;
        continue
    end

    y_plus = yy;
    y_plus(i) = yy(i) + t_plus / c(i);
    y_plus(j) = yy(j) - t_plus / c(j);

    y_minus = yy;
    y_minus(i) = yy(i) - t_minus / c(i);
    y_minus(j) = yy(j) + t_minus / c(j);

    % Convex weights: yy = theta*y_plus + (1-theta)*y_minus
    theta = t_minus / (t_plus + t_minus);

    w = W(k_split);

    % Replace point k_split with y_plus, add y_minus
    Y(:,k_split) = y_plus;
    W(k_split) = w*theta;

    Y = [Y, y_minus];
    W = [W; w*(1-theta)];
end

% Map each y-vertex to x-vertex
m = size(Y,2);
V = zeros(n,m);
for k = 1:m
    xx = l + d .* Y(:,k);
    % For fixed vars, d=0 so xx=l there automatically; but ensure exact:
    xx(~free) = l(~free);
    V(:,k) = xx;
end

lambda = W(:);

% Clean small weights and merge duplicates 
keep = lambda > 1e-15;
V = V(:,keep);
lambda = lambda(keep);

mergedV = [];
mergedL = [];
for k = 1:size(V,2)
    vk = V(:,k);
    idx = 0;
    for t = 1:size(mergedV,2)
        if norm(mergedV(:,t) - vk) <= 1e-12 * max(1,norm(vk))
            idx = t; break;
        end
    end
    if idx == 0
        mergedV = [mergedV, vk];
        mergedL = [mergedL; lambda(k)];
    else
        mergedL(idx) = mergedL(idx) + lambda(k);
    end
end
V = mergedV;
lambda = mergedL;

% Normalize
lambda(lambda<0 & lambda>-1e-14) = 0;
lambda = lambda / sum(lambda);

% Final check: x == V*lambda
xrec = V * lambda;
if norm(xrec - x) > 1e-8*max(1,norm(x))
    warning("Decomposition reconstruction error: %g", norm(xrec-x));
end

end

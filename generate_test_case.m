function generate_test_case(n, type)
% Create and RUN a test instance for Frank-Wolfe.
%
% Usage:
%   generate_test_case(n, "interior")
%   generate_test_case(n, "box_boundary")
%   generate_test_case(n, "active_linear")


fprintf("--------------------------------------------------\n");
fprintf(" Running test (%s), n = %d\n", type, n);
fprintf("--------------------------------------------------\n");

% boundary
l = zeros(n,1);
u = ones(n,1);

A = randn(n);
Q = A'*A + 0.1*eye(n);   % Q SPD
a = rand(n,1) + 0.1;     % positive

% Select type
switch lower(type)

    case 'interior'
        description = "Optimum inside the box, linear constraint inactive.";
        x_star = 0.3 + 0.4*rand(n,1);
        q = -2*Q*x_star;
        b = 0.1; 

    case 'box_boundary'
        description = "Optimum on the boundary of the box.";
        x_star = ones(n,1);   
        x_star(1) = 0.8;  
        q = -2*Q*x_star;
        b = 0.1;

    case 'active_linear'
        description = "Optimum lies on the linear constraint a'x = b.";
        x_star = rand(n,1);
        b = a' * x_star; % x_star is exactly on the linear constraint
        q = -2*Q*x_star;

    otherwise
        error("Test type not recognized.");
end

% Random starting point
x0 = rand(n,1);

% Set tolerance
eps = 1e-3;

% Run Frank Wolfe
[x_FW, f_FW, f_star] = frank_wolfe(Q, q, x0, a, b, l, u, eps);

% Print
fprintf("Description: %s\n", description);
fprintf("f(x_FW)     = %.6f\n", f_FW);
fprintf("f_star      = %.6f\n", f_star);
fprintf("--------------------------------------------------\n");

end

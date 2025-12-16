function run_tests(test_id)

% 2D tests only for visualization of the steps of the algorithm. 
% Two cases:
% 1) interior optimum
% 2) box boundary optimum
%
% Usage: run_tests(n), with n = 1,2

clc;
fprintf("Running test %d...\n", test_id);

switch test_id
    
    % ----------------------------------------------------
    case 1
        % TEST 1 — convergenza 
        % f(x) = x1^2 + x2^2 - 4 x1 - 2 x2
        Q = [1 0; 0 1];
        q = [-4; -2];

        a = [1; 1];
        b = 2;

        l = [0; 0];
        u = [3; 3];

        x0 = [2; 0]; % Admissible
        eps = 1e-3;

        [x_FW, f_FW, f_star] = frank_wolfe(Q, q, x0, a, b, l, u, eps);

        fprintf("x_FW = [%.6f ; %.6f]\n", x_FW(1), x_FW(2));
        fprintf("f(x) = %.6f\n", f_FW);
        fprintf("f_star = %.6f\n", f_star);

        return;
        
    % ----------------------------------------------------
    case 2
       % x_star is at boundary box, real minimum is (2,1) but it's not
       % feasible, best candidate is (1.9,1)
       Q = [1 0; 0 1];
       q = [-4; -2];

       a = [1; 1];
       b = 2;

       l = [0; 0];
       u = [1.9; 3];   

       x0 = [2; 0];         
       eps = 1e-3;

       [x_FW, f_FW, f_star] = frank_wolfe(Q, q, x0, a, b, l, u, eps);



       fprintf("x_FW = [%.6f ; %.6f]\n", x_FW(1), x_FW(2));
       fprintf("f(x_FW) = %.6f\n", f_FW);
       fprintf("f_star  = %.6f\n", f_star);


    % ----------------------------------------------------
    otherwise
        error("Test not valid");

       
end

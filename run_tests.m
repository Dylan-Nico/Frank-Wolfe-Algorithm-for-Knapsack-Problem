function run_tests(test_id)

clc;
fprintf("Running test %d...\n", test_id);

switch test_id
    
    % ----------------------------------------------------
    case 1
        % TEST 1 — convergenza lineare
        % f(x) = x1^2 + x2^2 - 4 x1 - 2 x2
        Q = [1 0; 0 1];
        q = [-4; -2];

        a = [1; 1];
        b = 2;

        l = [0; 0];
        u = [3; 3];

        x0 = [2; 0]; % ammissibile
        eps = 1e-6;

        [x_FW, f_FW, f_star] = frank_wolfe(Q, q, x0, a, b, l, u, eps);

        %fprintf("x_FW = "); disp(x_FW);
        %fprintf("a'*x = %.3f\n", a'*x_FW);
        fprintf("f(x) = %.6f\n", f_FW);
        fprintf("f_star = %.6f\n", f_star);

        return;
        
    % ----------------------------------------------------
    case 2
    n = 10;
    A = randn(n);
    Q = A'*A + 0.1*eye(n); % SPD
    q = randn(n,1)*0.1;

    l = zeros(n,1);
    u = ones(n,1);

    a = rand(n,1); 
    b = 2.5;

    x0 = rand(n,1);
    eps = 1e-6;

    [x_FW, f_FW, f_star] = frank_wolfe(Q, q, x0, a, b, l, u, eps);

    %fprintf("x_FW = "); disp(x_FW);
    %fprintf("a'*x = %.3f\n", a'*x_FW);
    fprintf("f(x) = %.6f\n", f_FW);
    fprintf("f_star = %.6f\n", f_star);


    % ----------------------------------------------------
    otherwise
        error("Test ID non valido (1-5).");

       
end

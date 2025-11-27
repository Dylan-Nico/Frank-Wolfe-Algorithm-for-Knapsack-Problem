% f(x) = x1^2 + x2^2 - 4*x1 - 2*x2
% dominio: 0 ≤ x ≤ 3
% vincolo:  x1 + x2 ≥ 2
% minimo noto: (2,1)
% Funzione convessa, partenza ammissibile

Q = [1 0; 0 1];
q = [-4; -2];

a = [1; 1];
b = 2;

l = [0; 0];
u = [3; 3];

% ---- PUNTO INIZIALE AMMISSIBILE ----
% a'*x0 = 2 >= b
x0 = [2; 0];

eps = 1e-8;

% ---- FRANK-WOLFE ----
[x_FW, f_FW] = frank_wolfe(Q, q, x0, a, b, l, u, eps);

fprintf("Soluzione FW:\n");
disp(x_FW);

fprintf("a'*x_FW = %.4f (>= 2)\n", a'*x_FW);
fprintf("f(x_FW) = %.6f\n", f_FW);

% Secondo test
% f(x) = x1^2 + x2^2 - 4*x1 - 2*x2
% dominio: 0 ≤ x ≤ 3
% vincolo:  x1 + x2 ≥ 2
% minimo noto: (2,1)


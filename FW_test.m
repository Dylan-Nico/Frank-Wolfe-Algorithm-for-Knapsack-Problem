clear; clc; close all;

% ----- Problem 1 -----
%Q = [4 -1; -1 2]; % Q well conditioned
%Q = [1e-5, 0; 0, 1]; % Q bad conditioned
Q = [4 -3; -3 2];
q = [-3; -1];

l = [0; 0];
u = [2; 2];

a = [1; 1];
b = 1;

% ----- Plot domain -----
figure; hold on; grid on;
plotBox2D_with_constraint(l, u, a, b);

x_uncon = -0.5 * (Q \ q);
plot(x_uncon(1), x_uncon(2), 'mo', 'MarkerSize',10,'MarkerFaceColor','m');
text(x_uncon(1)+0.05, x_uncon(2), 'Unconstrained minimum');


% ===== Plot delle curve di livello della funzione ===== %
[Xg, Yg] = meshgrid(linspace(l(1),u(1),200), linspace(l(2),u(2),200));

Fg = zeros(size(Xg));
for i = 1:numel(Xg)
    x_tmp = [Xg(i); Yg(i)];
    Fg(i) = x_tmp' * Q * x_tmp + q' * x_tmp;
end

contour(Xg, Yg, Fg, 20, 'LineColor',[0 0 0]);  % 20 curve – grigio semi-trasparente
    

title('Frank-Wolfe iterations (interactive)');

% ----- Start point (on boundary) -----
x = [0; 1];   % punto sul bordo
plot(x(1), x(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor','k');
text(x(1)+0.05, x(2), 'START');

disp('Premi INVIO per andare alla prossima iterazione');


for k = 1:30

    g = 2*Q*x + q;
    s = solveLP(g, a, b, l, u);

    d = s - x;

    denom = d'*(Q*d);
    if denom < 1e-12
        alpha = 1;
    else
        alpha = min(1, max(0, -(g'*d)/(2*denom)));
    end

    x_new = x + alpha*d;

    % --- Disegno della iterazione (solo dopo la prima!) ---
    plot([x(1), x_new(1)], [x(2), x_new(2)], 'k--');
    plot(x_new(1), x_new(2), 'ko', 'MarkerSize',5,'MarkerFaceColor','k');

    fprintf("Iter %d: x = [%.3f, %.3f]\n", k, x_new(1), x_new(2));

    pause; % <--- ASPETTA INVIO

    x = x_new;
end

% ----- Final point -----
plot(x(1), x(2), 'ro', 'MarkerSize',10,'MarkerFaceColor','r');
text(x(1)+0.05, x(2), 'OPTIMUM');

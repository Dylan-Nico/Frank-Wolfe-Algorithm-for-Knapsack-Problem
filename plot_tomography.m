function plot_tomography(Q, q, x, d, alpha)
% PLOT_TOMOGRAPHY Plotta la tomografia phi(gamma) = f(x + gamma d)
%   Q, q : definiscono f(x) = x'Qx + q'x
%   x    : punto corrente
%   d    : direzione FW (s - x)
%   alpha: optimum della line search esatta (punto rosso)

    % funzione obiettivo
    f = @(z) z' * Q * z + q' * z;

    % gamma discretizzati
    gammas = linspace(0, 1, 200);
    phi = zeros(size(gammas));

    % calcola phi(gamma)
    for i = 1:length(gammas)
        phi(i) = f(x + gammas(i)*d);
    end

    % plot
    figure; hold on; grid on;
    plot(gammas, phi, 'LineWidth', 1.8);
    xlabel('\gamma');
    ylabel('\phi(\gamma)');
    title('Tomografia della funzione lungo d');

    % punto della line search
    plot(alpha, f(x + alpha*d), 'ro', 'MarkerFaceColor', 'r');
    text(alpha, f(x + alpha*d), '  \alpha^*', 'Color','r','FontSize',12);
end

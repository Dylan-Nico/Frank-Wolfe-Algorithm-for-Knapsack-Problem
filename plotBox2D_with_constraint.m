function plotBox2D_with_constraint(l, u, a, b)

% Plots constraints of the problem in 2D cases (box + linear)

    % Create figure
    figure; hold on; grid on;

    % Plot box
    rectangle('Position', [l(1), l(2), u(1)-l(1), u(2)-l(2)], ...
              'EdgeColor', 'b', 'LineWidth', 2);

    % Plot the linear constraint
    x1 = linspace(l(1), u(1), 200);
    if abs(a(2)) > 1e-12
        x2 = (b - a(1)*x1) / a(2);
        plot(x1, x2, 'r', 'LineWidth', 2);
    else
        % vertical line
        x1_line = b / a(1);
        plot([x1_line x1_line], [l(2) u(2)], 'r', 'LineWidth', 2);
    end

    % Shade feasible region a^T x >= b
    [X,Y] = meshgrid(linspace(l(1),u(1),150), linspace(l(2),u(2),150));
    Z = a(1)*X + a(2)*Y >= b;

    % Shade
    contourf(X, Y, Z, [1 1], 'LineStyle', 'none');
    alpha(0.1);

    % labels
    xlabel('x_1');
    ylabel('x_2');
    title('Box constraint + linear constraint a^T x >= b');

    axis equal;
    xlim([l(1) u(1)]);
    ylim([l(2) u(2)]);

end

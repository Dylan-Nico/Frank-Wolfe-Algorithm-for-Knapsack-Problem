function FW_plot2D(Q,q,a,b,l,u,iterates)
    % iterates: matrice 2×K con gli x_k salvati

    figure; hold on; grid on;
    title('Frank-Wolfe iterations (2D)');
    xlabel('x_1'); ylabel('x_2');

    % ---- Plot dominio (box + vincolo) ----
    plotBox2D_with_constraint(l, u, a, b);

    % ---- Griglia per curve di livello ----
    [X, Y] = meshgrid(linspace(l(1),u(1),200), linspace(l(2),u(2),200));
    F = zeros(size(X));

    for i = 1:size(X,1)
        for j = 1:size(X,2)
            x = [X(i,j); Y(i,j)];
            F(i,j) = x'*Q*x + q'*x;
        end
    end

    % ---- Curve di livello ----
    contour(X, Y, F, 25, 'k'); 

    % ---- Iterazioni ----
    plot(iterates(1,:), iterates(2,:), 'o--', 'Color',[0 0 0], ...
        'MarkerFaceColor','k', 'LineWidth',1.4);

    % START
    text(iterates(1,1), iterates(2,1), ' START', 'FontWeight','bold');

    % END
    x_end = iterates(:,end);
    plot(x_end(1), x_end(2), 'ro', 'MarkerFaceColor','r');
    text(x_end(1), x_end(2)," OPTIMUM",'Color','r','FontWeight','bold');

end

function plot_gap(gaps)

    figure;
    semilogy(gaps, 'LineWidth', 1.8);
    xlabel('Iterazione');
    ylabel('Duality Gap (log scale)');
    title('Convergenza Frank–Wolfe');
    grid on;

end
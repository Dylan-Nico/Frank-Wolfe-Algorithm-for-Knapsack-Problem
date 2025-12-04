function plot_gaps(gaps, primal_errors)
%% Plots dual and primal gap

    % --- Figure ---
    figure;
    
    % Log scale plots
    semilogy(gaps, 'LineWidth', 1.5); hold on;
    semilogy(primal_errors, 'LineWidth', 1.5);
    
end

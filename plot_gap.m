function plot_gap2(gaps, primal_errors)

    k = min(length(gaps), length(primal_errors));
    gaps = gaps(1:k);
    primal_errors = primal_errors(1:k);
    % --- Messa in sicurezza ---
    gaps = gaps(:) + 1e-16;
    primal_errors = abs(primal_errors(:)) + 1e-16;
    k = length(gaps);

    % --- Soglie per gli zoom ---
    thresh_mid  = 1e-2;   % zoom centrale
    thresh_tail = 1e-3;   % zoom finale

    % Trovo gli indici da cui tagliare
    idx_mid  = find(primal_errors < thresh_mid, 1);
    if isempty(idx_mid), idx_mid = 1; end

    idx_tail = find(primal_errors < thresh_tail, 1);
    if isempty(idx_tail), idx_tail = floor(k*0.7); end

    % --- Figura ---
    figure;
    ylim([1e-4 1]);    
    
    % === (1) OVERVIEW COMPLETA ===
    subplot(3,1,1);
    semilogy(1:k, gaps, 'LineWidth', 1.5); hold on;
    semilogy(1:k, primal_errors, 'LineWidth', 1.5);
    grid on; title('Overview (All Iterations)');
    ylabel('Error');
    legend('Duality Gap', 'Primal Error');
    ax = gca;
    ax.Position = [0.12 0.15 0.8 0.75]; 
    
end

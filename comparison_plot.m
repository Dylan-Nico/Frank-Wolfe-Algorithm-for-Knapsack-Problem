function comparison_plot(gaps_FW, gaps_AW, pe_FW, pe_AW, type, n)
%% plots comparison between dual error of standard FW and Away FW (blue color)
%% and primal error of standard FW and Away FW (orange color)

% colors for mantaining consistency with singular plots
blue    = [0.0 0.2 0.6];   % dual gap
orange = [0.85 0.33 0.10]; % primal error


figure;
semilogy(1:length(gaps_FW), cummin(gaps_FW), 'LineWidth', 2.5, 'Color', blue, 'LineStyle', '--'); hold on;
semilogy(1:length(gaps_AW), cummin(gaps_AW), 'LineWidth', 2.5, 'Color', blue, 'LineStyle', '-');
grid on; xlabel('iteration'); ylabel('gap (your stopping gap)');
title(sprintf('Dual gap comparison (%s, n=%d)', type, n));
legend('FW standard','FW away');


figure;
semilogy(1:length(pe_FW), cummin(pe_FW), 'LineWidth', 2.5, 'Color', orange, 'LineStyle', '--'); hold on;
semilogy(1:length(pe_AW), cummin(pe_AW), 'LineWidth', 2.5, 'Color', orange, 'LineStyle', '-');
grid on; xlabel('iteration'); ylabel('primal error');
title(sprintf('Primal error comparison (%s, n=%d)', type, n));
legend('FW standard','FW away');

% Data (MPa) vs Temperature (K)
T  = [293 193  93];    % K (given, descending order is OK for pchip/makima)
Xt = [3298 3183 3168]; % MPa
Xc = [1753 1732 1715]; % MPa
Yt = [  50   50   53]; % MPa
Yc = [  96  148  259]; % MPa
S12= [  81  126  169]; % MPa

% Query temperatures for plotting
Tq = linspace(20, 293, 400); % includes extrapolation down to 20 K
T20 = 20;

% Helper function for plotting one property
plot_compare = @(name, y) plot_one(name, T, y, Tq, T20);

% Plot each property in its own figure (clear & readable)
plot_compare('X_t (MPa)', Xt);
plot_compare('X_c (MPa)', Xc);
plot_compare('Y_t (MPa)', Yt);
plot_compare('Y_c (MPa)', Yc);
plot_compare('S_{12} (MPa)', S12);

% ---------- Local function ----------
function plot_one(ylabelStr, T, y, Tq, T20)

    % Interpolate / extrapolate
    y_pchip  = pchip(T, y, Tq);
    y_makima = makima(T, y, Tq);

    % Values at 20 K
    y20_pchip  = pchip(T, y, T20);
    y20_makima = makima(T, y, T20);

    % Plot
    figure; hold on; grid on;
    plot(Tq, y_pchip,  '-',  'LineWidth', 1.8, 'DisplayName', 'pchip');
    plot(Tq, y_makima, '--', 'LineWidth', 1.8, 'DisplayName', 'makima');

    % Original data points
    plot(T, y, 'o', 'MarkerSize', 7, 'LineWidth', 1.5, 'DisplayName', 'data');

    % Mark 20 K estimates
    plot(T20, y20_pchip,  's', 'MarkerSize', 8, 'LineWidth', 1.6, ...
        'DisplayName', sprintf('pchip @20K = %.2f', y20_pchip));
    plot(T20, y20_makima, 'd', 'MarkerSize', 8, 'LineWidth', 1.6, ...
        'DisplayName', sprintf('makima @20K = %.2f', y20_makima));

    xlabel('Temperature (K)');
    ylabel(ylabelStr);
    title(['Interpolation / extrapolation comparison: ', ylabelStr]);

    % Nice legend placement
    legend('Location','best');

    % Print to console for quick comparison
    fprintf('%s:\n', ylabelStr);
    fprintf('  pchip  @20K = %.4f\n', y20_pchip);
    fprintf('  makima @20K = %.4f\n\n', y20_makima);

end

function plot_pareto_front(csvfile)

    % -------------------------------
    % Load CSV
    % -------------------------------
    T = readtable(csvfile);

    % Column names (exact as given)
    pack = T.("PackingEff");      % Packing efficiency
    stress = T.max_vm;            % Von Mises stress

    % -------------------------------
    % Sphere and Cylinder reference values
    % -------------------------------
    pack_sphere = 0.5236123585;
    stress_sphere = 54.26379776;

    pack_cyl = 0.69;
    stress_cyl = 69.0;

    % -------------------------------
    % Plot all points
    % -------------------------------
    figure; hold on;
    scatter(pack, stress, 55, 'filled', ...
        'MarkerFaceColor',[0.6 0.6 0.9], 'MarkerEdgeColor','k');

    % Highlight sphere
    scatter(pack_sphere, stress_sphere, 120, 'g', 'filled', ...
        'MarkerEdgeColor','k', 'LineWidth',1.2);
    % text(pack_sphere+0.005, stress_sphere, 'Sphere', 'FontSize',15, 'Color','g');

    % Highlight cylinder
    scatter(pack_cyl, stress_cyl, 120, 'm', 'filled', ...
        'MarkerEdgeColor','k', 'LineWidth',1.2);
    % text(pack_cyl+0.005, stress_cyl, 'Cylinder', 'FontSize',15, 'Color','m');

    xlabel('Packing Efficiency');
    ylabel('Max von Mises Stress (MPa)');
    title('Pareto Front: Stress vs Packing Efficiency');
    grid on;

    % -------------------------------
    % Compute Pareto front
    % -------------------------------
    data = [pack stress];
    isPareto = find_pareto(data);

    % Extract Pareto points
    pack_p = pack(isPareto);
    stress_p = stress(isPareto);

    % Sort for smooth line
    [pack_p, I] = sort(pack_p);
    stress_p = stress_p(I);

    % -------------------------------
    % Plot Pareto curve
    % -------------------------------
    plot(pack_p, stress_p, 'r-', 'LineWidth', 2);
    scatter(pack_p, stress_p, 75, 'r', 'filled');

    legend('All designs', 'Sphere', 'Cylinder', 'Pareto optimal', ...
        'Location','best');

    hold off;

end


% -------------------------------------------------------
% Pareto helper function
% -------------------------------------------------------
function isPareto = find_pareto(data)
    N = size(data,1);
    isPareto = true(N,1);

    for i = 1:N
        for j = 1:N
            if j == i, continue; end

            better_or_equal = ...
                data(j,1) >= data(i,1) && ...   % packing higher or equal
                data(j,2) <= data(i,2);         % stress lower or equal

            strictly_better = ...
                data(j,1) > data(i,1) || ...
                data(j,2) < data(i,2);

            if better_or_equal && strictly_better
                isPareto(i) = false;
                break;
            end
        end
    end
end

% Run
plot_pareto_front("cleaned_file.csv");

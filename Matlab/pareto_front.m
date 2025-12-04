function plot_pareto_front(csvfile)

    % -------------------------------
    % Load CSV
    % -------------------------------
    T = readtable(csvfile);


    % Column names (exact as given)
    pack = T.("PackingEff");  % Packing efficiency
    stress = T.max_vm;         % Von Mises stress

    % -------------------------------
    % Plot all points
    % -------------------------------
    figure; hold on;
    scatter(pack, stress, 55, 'filled', ...
        'MarkerFaceColor',[0.6 0.6 0.9], 'MarkerEdgeColor','k');

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

    legend('All designs', 'Pareto optimal', 'Location','best');
    hold off;

end



% -------------------------------------------------------
% Pareto helper function
% -------------------------------------------------------
function isPareto = find_pareto(data)
    % data(:,1) = packing efficiency (maximize)
    % data(:,2) = stress (minimize)
    %
    % A point i is Pareto-optimal if:
    %   No point j has:
    %      pack_j >= pack_i AND stress_j <= stress_i
    %      and at least one strict improvement

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


% plot_pareto_front("combined_pareto.csv");
plot_pareto_front("cleaned_file.csv");

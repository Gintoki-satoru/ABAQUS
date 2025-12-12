function plot_stress_heatmap(csvfile)

    T = readtable(csvfile);

    % Extract columns
    n1_vals = T.n1;
    n2_vals = T.n2;
    stress  = T.max_vm;

    % Unique sorted values for combined pareto
    % n1_unique = unique(n1_vals);
    % n2_unique = unique(n2_vals);
    % -----------------------------------------------
    % FORCE n1 and n2 to be this exact range:
    % 1.0, 0.95, 0.90, ..., 0.60
    % -----------------------------------------------
    n_unique = 1.0:-0.05:0.6;
    n1_unique = n_unique;
    n2_unique = n_unique;

    % Initialize heatmap matrix
    H = nan(length(n2_unique), length(n1_unique));

    % Fill matrix (choose MIN stress for duplicate (n1,n2))
    for i = 1:length(n2_unique)
        for j = 1:length(n1_unique)
            
            idx = find(n1_vals == n1_unique(j) & n2_vals == n2_unique(i));

            if ~isempty(idx)
                % use min stress across all matching geometries
                H(i,j) = min(stress(idx));
            end
        end
    end

    % Plot heatmap
    figure;
    imagesc(n1_unique, n2_unique, H);
    colorbar;
    xlabel('n1');
    ylabel('n2');
    title('Stress Heatmap: min stress for each (n1,n2)');
    set(gca, 'YDir', 'normal');

    % Add numbers on heatmap
    for i = 1:size(H,1)
        for j = 1:size(H,2)
            if ~isnan(H(i,j))
                text(n1_unique(j), n2_unique(i), ...
                     sprintf('%.1f', H(i,j)), ...
                     'HorizontalAlignment','center', ...
                     'Color','white','FontSize',10);
            end
        end
    end
end

% plot_stress_heatmap("combined_pareto.csv")
plot_stress_heatmap("cleaned_file.csv")

%% heat map for packing efficiency

function plot_packing_heatmap(csvfile)

    T = readtable(csvfile);

    % Extract columns
    n1_vals = T.n1;
    n2_vals = T.n2;
    pack    = T.("PackingEff") ;

    % -----------------------------------------------
    % Force n-grid: 1.0, 0.95, ..., 0.60
    % -----------------------------------------------
    n_unique = 1.0:-0.05:0.6;
    n1_unique = n_unique;
    n2_unique = n_unique;

    % Heatmap matrix (packing efficiency)
    H = nan(length(n2_unique), length(n1_unique));

    % Fill matrix using MEAN packing efficiency for duplicates
    for i = 1:length(n2_unique)
        for j = 1:length(n1_unique)

            idx = find(abs(n1_vals - n1_unique(j)) < 1e-6 & ...
                       abs(n2_vals - n2_unique(i)) < 1e-6);

            if ~isempty(idx)
                H(i,j) = max(pack(idx));   % mean PE across geometry variations
            end
        end
    end

    % -----------------------------------------------
    % Plot heatmap
    % -----------------------------------------------
    figure;
    imagesc(n1_unique, n2_unique, H);
    colormap("turbo");  % colorful & clear
    colorbar;

    xlabel('n1');
    ylabel('n2');
    title('Packing Efficiency Heatmap over (n1, n2)');
    set(gca, 'YDir', 'normal', 'FontSize', 12);

    % Add labels (formatted %)
    for i = 1:size(H,1)
        for j = 1:size(H,2)
            if ~isnan(H(i,j))
                text(n1_unique(j), n2_unique(i), ...
                     sprintf('%.2f', H(i,j)), ...
                     'HorizontalAlignment','center', ...
                     'Color','white','FontSize',10);
            end
        end
    end
end

% Call the function
plot_packing_heatmap("cleaned_file.csv");

%% compare packing efficiency

function compare_packing_n1_n2(csvfile, n_fixed)

    if nargin < 2
        n_fixed = 0.8;   % default representative value
    end

    T = readtable(csvfile);

    % --- Extract relevant columns ---
    n1_vals = T.n1;
    n2_vals = T.n2;

    % auto-detect packing efficiency column
    col_pack = find(contains(T.Properties.VariableNames, "pack", "IgnoreCase", true), 1);
    if isempty(col_pack)
        error("Could not find packing efficiency column (name containing 'pack').");
    end
    pack = T{:, col_pack};

    % --- Define exponent grid: 1.0, 0.95, ..., 0.60 ---
    n_unique = 1.0:-0.05:0.6;
    n1_unique = n_unique;
    n2_unique = n_unique;

    % --- Build packing efficiency matrix H_pack(i,j) ---
    H_pack = nan(length(n2_unique), length(n1_unique));

    for i = 1:length(n2_unique)
        for j = 1:length(n1_unique)
            idx = find( abs(n1_vals - n1_unique(j)) < 1e-6 & ...
                        abs(n2_vals - n2_unique(i)) < 1e-6 );
            if ~isempty(idx)
                H_pack(i,j) = mean(pack(idx));  % mean packing efficiency for this (n1,n2)
            end
        end
    end

    % --- Find nearest index to n_fixed on grid ---
    [~, idx_fix_n1] = min(abs(n1_unique - n_fixed));  % column index (n1 = n_fixed)
    [~, idx_fix_n2] = min(abs(n2_unique - n_fixed));  % row index (n2 = n_fixed)

    n_fixed_n1 = n1_unique(idx_fix_n1);
    n_fixed_n2 = n2_unique(idx_fix_n2);

    % --- Curve 1: n2 fixed = n_fixed, vary n1 ---
    curve_n1 = H_pack(idx_fix_n2, :);   % row: fixed n2
    % --- Curve 2: n1 fixed = n_fixed, vary n2 ---
    curve_n2 = H_pack(:, idx_fix_n1);   % column: fixed n1

    % --- Plot both on same axes ---
    figure; hold on;
    plot(n1_unique, curve_n1, 'o-b', 'LineWidth', 2, 'DisplayName', ...
         sprintf('n2 = %.2f, n1 varying', n_fixed_n2));
    plot(n2_unique, curve_n2, 's-r', 'LineWidth', 2, 'DisplayName', ...
         sprintf('n1 = %.2f, n2 varying', n_fixed_n1));

    xlabel('Exponent value');
    ylabel('Packing Efficiency');
    title(sprintf('Effect of decreasing n1 vs n2 on Packing Efficiency (n* = %.2f)', n_fixed));
    grid on;
    legend('Location','best');
    set(gca, 'XDir','reverse');  % so 1.0 → 0.6 goes left to right if you like

    hold off;

    % --- Simple quantitative comparison: average changes ---
    % Drop NaNs so we only use valid points
    valid1 = ~isnan(curve_n1);
    valid2 = ~isnan(curve_n2);

    % sort by exponent (from 1 → 0.6) for consistent Δ
    [n_sorted1, i1] = sort(n1_unique(valid1), 'descend');
    c1 = curve_n1(valid1);
    c1 = c1(i1);

    [n_sorted2, i2] = sort(n2_unique(valid2), 'descend');
    c2 = curve_n2(valid2);
    c2 = c2(i2);

    if numel(c1) > 1 && numel(c2) > 1
        dpack1 = c1(end) - c1(1);   % change in packing when n1 drops from ~1 → ~0.6
        dpack2 = c2(end) - c2(1);   % change in packing when n2 drops from ~1 → ~0.6

        fprintf('\nChange in packing eff when n1 drops (n2 fixed = %.2f):  %.4f\n', ...
                n_fixed_n2, dpack1);
        fprintf('Change in packing eff when n2 drops (n1 fixed = %.2f):  %.4f\n\n', ...
                n_fixed_n1, dpack2);
    end

end

compare_packing_n1_n2("cleaned_file.csv", 0.8);   % or 1.0, 0.9, etc.

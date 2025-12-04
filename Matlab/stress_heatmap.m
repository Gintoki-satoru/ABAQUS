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
                H(i,j) = mean(pack(idx));   % mean PE across geometry variations
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

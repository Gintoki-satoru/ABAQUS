function plot_slenderness(csvfile)

    % ------------------------------------------
    % Load the table
    % ------------------------------------------
    T = readtable(csvfile);

    % Auto-detect columns
    col_a = find_col(T, "a_mm");
    col_b = find_col(T, "b_mm");
    col_c = find_col(T, "c_mm");
    col_stress = find_col(T, "max_vm");

    a = double(T{:,col_a});
    b = double(T{:,col_b});
    c = double(T{:,col_c});
    stress = double(T{:,col_stress});

    % ------------------------------------------
    % PLOT A: Axisymmetric family (a ≈ b)
    % ------------------------------------------
    idx_axial = abs(a - b) < 1e-6;   % almost equal

    lambda1 = c(idx_axial) ./ a(idx_axial);
    stress1 = stress(idx_axial);

    figure;
    scatter(lambda1, stress1, 60, 'filled');
    xlabel('\lambda_1 = c / a   (axisymmetric: a = b)');
    ylabel('Stress (MPa)');
    title('Stress vs Slenderness Ratio c/a (for a = b)');
    grid on;


    % ------------------------------------------
    % PLOT B: Sideways family (a ≈ c)
    % ------------------------------------------
    idx_side = abs(a - c) < 1e-6;

    lambda2 = b(idx_side) ./ a(idx_side);
    stress2 = stress(idx_side);

    figure;
    scatter(lambda2, stress2, 60, 'filled');
    xlabel('\lambda_2 = b / a   (sideways: a = c)');
    ylabel('Stress (MPa)');
    title('Stress vs Lateral Ratio b/a (for a = c)');
    grid on;

end


% ------------------------------------------
% Helper: find a column from partial key
% ------------------------------------------
function idx = find_col(T, key)
    idx = find(contains(T.Properties.VariableNames, key, 'IgnoreCase', true), 1);
    if isempty(idx)
        error("Column containing '%s' not found. Columns are:\n%s\n", ...
              key, strjoin(T.Properties.VariableNames, ", "));
    end
end
plot_slenderness("combined_pareto.csv");

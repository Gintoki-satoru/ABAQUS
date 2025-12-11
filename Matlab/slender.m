function analyze_triaxial_geometries(csvfile)

    % ---- Load data ----
    T = readtable(csvfile);

    a = T.a_mm;
    b = T.b_mm;
    c = T.c_mm;
    n1 = T.n1;
    n2 = T.n2;
    stress = T.max_vm;

    tol = 1e-3;   % tolerance for saying two axes are "equal"

    % ---- Classification: symmetric / axisymmetric / triaxial ----
    eq_ab = abs(a - b) < tol;
    eq_ac = abs(a - c) < tol;
    eq_bc = abs(b - c) < tol;

    is_symmetric   = eq_ab & eq_ac;                     % a ≈ b ≈ c
    is_axisym_ab   = eq_ab & ~eq_ac & ~eq_bc;           % a ≈ b ≠ c
    is_axisym_ac   = eq_ac & ~eq_ab & ~eq_bc;           % a ≈ c ≠ b
    is_axisym_bc   = eq_bc & ~eq_ab & ~eq_ac;           % b ≈ c ≠ a

    is_triaxial    = ~(is_symmetric | is_axisym_ab | is_axisym_ac | is_axisym_bc);

    % Only two groups for now: triaxial vs others
    group = strings(height(T),1);
    group(is_triaxial)      = "Unsymmetric (a≠b≠c)";
    group(~is_triaxial)     = "symmetric / axisymmetric (a=b)";

    % ---- 1) Boxplot: stress distribution per group ----
    figure;
    boxplot(stress, group);
    ylabel('Max von Mises stress (MPa)');
    title('Stress comparison: Unsymmetric vs symmetric/axisymmetric geometries');
    grid on;

    % ---- 2) Check per (n1,n2) whether triaxial ever beats others ----
    n1_unique = unique(n1);
    n2_unique = unique(n2);

    bad_pairs = 0;
    for i = 1:length(n2_unique)
        for j = 1:length(n1_unique)

            mask = abs(n1 - n1_unique(j)) < 1e-6 & ...
                   abs(n2 - n2_unique(i)) < 1e-6;

            if ~any(mask), continue; end

            tri_mask   = mask & is_triaxial;
            other_mask = mask & ~is_triaxial;

            % only compare if both families are present for this (n1,n2)
            if any(tri_mask) && any(other_mask)
                tri_min   = min(stress(tri_mask));
                other_min = min(stress(other_mask));

                if tri_min <= other_min
                    bad_pairs = bad_pairs + 1;
                    fprintf('Warning: at (n1=%.3f, n2=%.3f), tri_min <= other_min (%.2f <= %.2f)\n', ...
                            n1_unique(j), n2_unique(i), tri_min, other_min);
                end
            end
        end
    end

    if bad_pairs == 0
        fprintf('\nResult: For every (n1,n2) where both types exist, the BEST triaxial shape has HIGHER stress than the BEST symmetric/axisymmetric shape.\n');
    else
        fprintf('\nResult: %d exponent pairs where triaxial was not worse – check above warnings.\n', bad_pairs);
    end
end

analyze_triaxial_geometries("cleaned_file.csv");

%%

function plot_stress_vs_slenderness(csvfile)

    % Load data
    T = readtable(csvfile);

    % --- Restrict to n-range 0.6 to 1.0 ---
    T = T(T.n1 >= 0.6 & T.n1 <= 1.0 & ...
          T.n2 >= 0.6 & T.n2 <= 1.0, :);

    % Extract
    a = T.a_mm;
    b = T.b_mm;
    c = T.c_mm;

    % Only keep axisymmetric shapes (a = b)
    mask_ab = abs(a - b) < 1e-6;
    T = T(mask_ab, :);

    % Compute slenderness λ = c/a
    lambda = T.c_mm ./ T.a_mm;

    % Remove undesired λ < 1 values
    T = T(lambda >= 1, :);
    lambda = lambda(lambda >= 1);

    % ----------------------------------------------------
    % PLOT 1: n2 = 1, n1 varies, λ = [1,3,5,7,9]
    % ----------------------------------------------------

    slender_target_1 = [1 3 5 7 9];
    n2_fixed = 1;
    n1_values = unique(T.n1);

    figure; hold on;

    for n1_val = n1_values'

        mask = abs(T.n1 - n1_val) < 1e-6 & abs(T.n2 - n2_fixed) < 1e-6;

        if ~any(mask), continue; end

        lam = lambda(mask);
        s   = T.max_vm(mask);

        % Select exact slenderness targets (rounded)
        keep = ismember(round(lam), slender_target_1);
        lam  = lam(keep);
        s    = s(keep);

        if numel(lam) < 2, continue; end

        % Sort before interpolation
        [lam_sorted, idx] = sort(lam);
        s_sorted = s(idx);

        % Spline interpolation for smooth curve
        lam_fine = linspace(min(lam_sorted), max(lam_sorted), 200);
        s_fine = spline(lam_sorted, s_sorted, lam_fine);

        plot(lam_fine, s_fine, 'LineWidth', 1.5);
        plot(lam_sorted, s_sorted, 'o', 'DisplayName', sprintf('n1 = %.2f', n1_val));
    end

    xlabel('Slenderness ratio \lambda = c/a')
    ylabel('Stress (MPa)')
    title('Stress vs Slenderness ratio (n2 = 1, n1 decreasing)')
    legend show; grid on;
    hold off;


    % ----------------------------------------------------
    % PLOT 2: n1 = 1, n2 varies, λ = [1,3,5]
    % ----------------------------------------------------

    slender_target_2 = [1 3 5];
    n1_fixed = 1;
    n2_values = unique(T.n2);

    figure; hold on;

    for n2_val = n2_values'

        mask = abs(T.n1 - n1_fixed) < 1e-6 & abs(T.n2 - n2_val) < 1e-6;

        if ~any(mask), continue; end

        lam = lambda(mask);
        s   = T.max_vm(mask);

        keep = ismember(round(lam), slender_target_2);
        lam = lam(keep); s = s(keep);

        if numel(lam) < 2, continue; end

        % Sort for spline
        [lam_sorted, idx] = sort(lam);
        s_sorted = s(idx);

        % Spline interpolation
        lam_fine = linspace(min(lam_sorted), max(lam_sorted), 200);
        s_fine = spline(lam_sorted, s_sorted, lam_fine);

        plot(lam_fine, s_fine, 'LineWidth', 1.5);
        plot(lam_sorted, s_sorted, 's', 'DisplayName', sprintf('n2 = %.2f', n2_val));
    end

    xlabel('Slenderness ratio \lambda = c/a')
    ylabel('Stress (MPa)')
    title('Stress vs Slenderness ratio (n1 = 1, n2 decreasing)')
    legend show; grid on;
    hold off;

end



plot_stress_vs_slenderness("cleaned_file.csv");

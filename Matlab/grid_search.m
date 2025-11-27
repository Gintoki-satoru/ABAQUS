function results = generate_superellipsoid_abc_n()

    % ---- INPUTS ----
    n_vals = [1.0, 0.8, 0.6, 0.4];
    ratio_vals = [1, 3, 5];

    V_inner_target = 0.0587e9;      
    V_wall_target  = 0.0017496e9;   
    V_outer_target = V_inner_target + V_wall_target;

    % ---- Generate unique (a,b) by enforcing a <= b ----
    ab_pairs = [];
    for ai = 1:numel(ratio_vals)
        for bi = 1:numel(ratio_vals)
            a = ratio_vals(ai);
            b = ratio_vals(bi);
            if a <= b      % enforce uniqueness of (a,b)
                ab_pairs = [ab_pairs; a, b];
            end
        end
    end
    
    % How many?
    % disp(ab_pairs)  % should show 6 unique rows

    % ---- All c-values ----
    c_vals = ratio_vals(:);  % 1, 3, 5

    % Total geometric combinations = 6 ab-pairs * 3 c-values = 18
    % Total with n = 18 * 5 = 90

    total_rows = size(ab_pairs,1) * numel(c_vals) * numel(n_vals);
    data = zeros(total_rows, 9); 
    % columns: n, a_ratio, b_ratio, c_ratio, a, b, c, t, V_out

    idx = 1;

    % ---- MAIN LOOP ----
    for ni = 1:numel(n_vals)
        n = n_vals(ni);

        for abi = 1:size(ab_pairs,1)
            a_ratio = ab_pairs(abi,1);
            b_ratio = ab_pairs(abi,2);

            for ci = 1:numel(c_vals)
                c_ratio = c_vals(ci);

                % --- Step 1: base shape before scaling ---
                a0 = a_ratio;
                b0 = b_ratio;
                c0 = c_ratio;

                % Compute base volume
                V0 = superellipsoid_volume(a0, b0, c0, n, n);

                % Scale factor to match inner volume
                s = (V_inner_target / V0)^(1/3);

                % Scaled inner dimensions
                a_in = s*a0;
                b_in = s*b0;
                c_in = s*c0;

                % --- Step 2: solve thickness t so outer volume matches ---
                fun = @(t) superellipsoid_volume(a_in+t, b_in+t, c_in+t, n, n) - V_outer_target;
                t_guess = 5; 
                options = optimset('Display','off');
                t = fzero(fun, [0, 5], options);

                V_out = superellipsoid_volume(a_in+t, b_in+t, c_in+t, n, n);

                % --- Store ---
                data(idx,:) = [n, a_ratio, b_ratio, c_ratio, a_in, b_in, c_in, t, V_out];
                idx = idx + 1;
            end
        end
    end

    % ---- OUTPUT TABLE ----
    results = array2table(data, ...
         'VariableNames', ...
         {'n','ratio_a','ratio_b','ratio_c','a','b','c','t','V_outer'});
    writetable(results, 'results.csv');
end



% ---- Your volume function ----
function V = superellipsoid_volume(a, b, c, n1, n2)
    e1 = 2/n1;
    e2 = 2/n2;
    V = 8 * a * b * c * ...
        (gamma(1 + 1/e1)^2 * gamma(1 + 1/e2)) / ...
        (gamma(1 + 2/e1) * gamma(1 + (1/e2 + 2/e1)));
end

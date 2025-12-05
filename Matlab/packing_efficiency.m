function V = superellipsoid_volume(a, b, c, n1, n2)
    e1 =2/n1;
    e2 =2/n2;
    V = 8 * a * b * c * (gamma(1 + 1/e1)^2 * gamma(1 + 1/e2)) / ...
        (gamma(1 + 2/e1) * gamma(1 + (1/e2 + 2/e1)));
end

a = 184;
b = 184;
c = 552;
n1 = 1;
n2 = 2;
thick = 2.33;

V_inner = superellipsoid_volume(a, b, c, n1, n2);

V_cube = 8*a*b*c;

eff = 5.870308322012401e+07/V_cube*100;
fprintf('Packing Efficiency = %.3f \n', eff);

%%
function plot_packing_n1_n2_range()

    % -----------------------------
    % Define exponent range
    % -----------------------------
    n_range = 0.4 : 0.1 : 2;       % from 1 to 2
    n1_vals = n_range;
    n2_vals = n_range;

    % -----------------------------
    % Initialize packing efficiency matrix
    % -----------------------------
    H = nan(length(n2_vals), length(n1_vals));

    % -----------------------------
    % Loop through n1, n2 pairs
    % -----------------------------
    for i = 1:length(n2_vals)
        for j = 1:length(n1_vals)

            n1 = n1_vals(j);
            n2 = n2_vals(i);

            V = superellipsoid_volume(1, 1, 1, n1, n2);
            H(i,j) = V / 8;   % packing efficiency
        end
    end

    % -----------------------------
    % Plot heatmap
    % -----------------------------
    figure;
    imagesc(n1_vals, n2_vals, H);
    set(gca, 'YDir','normal');
    colormap turbo; colorbar;

    xlabel('n1'); ylabel('n2');
    title('Packing Efficiency for n1,n2 = 1 to 2 (a=b=c=1)');
    
    % Add numeric labels
    for i = 1:size(H,1)
        for j = 1:size(H,2)
            text(n1_vals(j), n2_vals(i), sprintf('%.2f', H(i,j)), ...
                'Color','white', 'FontSize',8, ...
                'HorizontalAlignment','center');
        end
    end
end


% -----------------------------
% Superellipsoid Volume Function
% -----------------------------

plot_packing_n1_n2_range();

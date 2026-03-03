function V = superellipsoid_volume(a1, a2, a3, n1, n2)
    eps1 = n1;   % controls z-related rounding
    eps2 = n2;   % controls xy-plane rounding

    % Beta function (MATLAB has beta(); otherwise use gamma relation)
    B = @(p,q) beta(p,q);
    % If you prefer no beta(), uncomment and use:
    % B = @(p,q) gamma(p).*gamma(q)./gamma(p+q);

    V = 2*a1*a2*a3 * eps1*eps2 * ...
        B(eps1/2 + 1, eps1) * ...
        B(eps2/2, eps2/2);
end
a = 46.05;
b = 46.05;
c = 120*a;
n1 = 0.45;
n2 = 1;
thick = 1.503;

V_inner = superellipsoid_volume(a, b, c, n1, n2);

V_cube = 8*a*b*c;

eff = V_inner/V_cube*100;
fprintf('Packing Efficiency = %.3f \n', eff);

%%
function plot_packing_n1_n2_range()

    % -----------------------------
    % Define exponent range
    % -----------------------------
    n_range = 0.1 : 0.1 : 1;       % from 1 to 2
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

    xlabel('n_1'); ylabel('n_2');
    title('Packing Efficiency for n_1,n_2 = 0.1 to 1 (a=b=c=1)');
    
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

%%

% Geometry (can be arbitrary since efficiency is shape-only)
a = 241.09;
b = 241.09;
c = 241.09;

% Fixed exponent
n2 = 1;

% Vary n1
n1_vals = linspace(1, 0.1, 100);

packing_eff = zeros(size(n1_vals));

% Bounding box volume
V_box = 8 * a * b * c;

% Loop over n1
for i = 1:length(n1_vals)
    n1 = n1_vals(i);
    V_se = superellipsoid_volume(a, b, c, n1, n2);
    packing_eff(i) = V_se / V_box;
end

% ---- Reference geometries in the same bounding box ----
rs = 241.09;                 % radius limited by half-width (a=b assumed)
              % total available height in box
V_box = 8 *rs^3;
% Sphere (fits when c>=a; here c=a)
V_sphere = (4/3) * pi * rs^3;
PE_sphere = V_sphere / V_box;

r = 152.5;
H = 905; 
V_box = H*4*r^2;
% Cylinder + hemispherical domes ("capsule") with total height H
L = max(0, H - 2*r);   % cylinder length (>=0)
V_capsule = pi*r^2*L + (4/3)*pi*r^3;
PE_capsule = V_capsule / V_box;

% Plot
figure;
plot(n1_vals, packing_eff, 'LineWidth', 2); hold on;
% yline(PE_sphere, '--', 'Spherical geometry', 'LineWidth', 1.5);
% yline(PE_capsule, ':', 'Cylindrical geometry', 'LineWidth', 1.5);
hold off;

xlabel('n_1');
ylabel('Packing Efficiency');
title('Packing Efficiency vs n_1 (n_2 = 1)');
grid on;
% legend('Superellipsoid', 'Sphere', 'Cyl + hemis domes', 'Location', 'best');

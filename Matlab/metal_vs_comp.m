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
c = 5526.54;
n1 = 0.45;
n2 = 1;
thick_comp = 0.64;

V_inner_comp = superellipsoid_volume(a, b, c, n1, n2);

V_outer_comp = superellipsoid_volume(a+thick_comp, b+thick_comp, c+thick_comp, n1, n2);

V_material_comp = V_outer_comp - V_inner_comp;

rho_f = 1810;     % fiber density (kg/m^3)
rho_m = 1310;     % matrix density (kg/m^3)
Vf    = 0.6;      % fiber volume fraction

% ---- Composite density ----
rho_c = Vf * rho_f + (1 - Vf) * rho_m;

fprintf('Composite density = %.2f kg/m^3\n', rho_c);

thick_m = 0.509;

V_inner_m = superellipsoid_volume(a, b, c, n1, n2);

V_outer_m = superellipsoid_volume(a+thick_m, b+thick_m, c+thick_m, n1, n2);

V_material_m = V_outer_m - V_inner_m;

rho_m = 2820;

% ---- Mass values (replace with your actual values) ----
mass_composite = V_material_comp* 1e-9*rho_c;
mass_metal = V_material_m* 1e-9*rho_m;

data = [mass_metal mass_composite];

figure
b = bar(data, 0.5);

grid on
ylabel('Tank Mass (kg)')
title('Mass Comparison for Super-Ellipsoid (n_1 = 0.45) Geometry')

set(gca, 'XTickLabel', {'Metal Tank','Composite Tank'})

% ---- Percentage weight saving ----
weight_saving = (mass_metal - mass_composite)/mass_metal * 100;

text(2, mass_composite, ...
     sprintf('%.1f%% lighter', weight_saving), ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', ...
     'FontSize',10);
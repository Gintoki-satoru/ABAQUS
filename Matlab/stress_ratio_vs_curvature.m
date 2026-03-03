Sphere = readmatrix('hoopvmeri.xlsx', 'Range', 'B:E');
super_15 = readmatrix('hoopvmeri.xlsx', 'Range', 'G:J');
super_75 = readmatrix('hoopvmeri.xlsx', 'Range', 'L:O');

% ---- Stress ratio from FE ----
sphere_ratio = super_75(:,4) ./ super_75(:,3);
z_data = super_75(:,2);

% ---- Super-ellipse parameters ----
a = 103.59;    % radial semi-axis (mm)
c = 932.3;    % longitudinal semi-axis (mm)
n = 0.65;      % exponent

theta = linspace(0, pi/2, 2000);

% ---- Geometry ----
r_geom = a .* (cos(theta).^n);
z_geom = c .* (sin(theta).^n);

% ---- First derivatives ----
dr_dth = -a * n .* (cos(theta).^(n-1)) .* sin(theta);
dz_dth =  c * n .* (sin(theta).^(n-1)) .* cos(theta);

% ---- Second derivatives ----
d2r_dth2 = -a*n*((n-1)*(cos(theta).^(n-2)).*(-sin(theta).^2) ...
                 + (cos(theta).^(n-1)).*cos(theta));

d2z_dth2 =  c*n*((n-1)*(sin(theta).^(n-2)).*(cos(theta).^2) ...
                 - (sin(theta).^(n-1)).*sin(theta));

% ---- Curvature ----
kappa = abs(dr_dth .* d2z_dth2 - dz_dth .* d2r_dth2) ./ ...
        ( (dr_dth.^2 + dz_dth.^2).^(3/2) );

% ---- Scale geometry to match FE axis ----
z_geom_scaled = z_geom / max(z_geom) * max(z_data);

% ---- Normalize curvature for visualization ----
kappa_scaled = kappa / max(kappa) * max(r_geom);

figure

% ---- LEFT: stress ratio ----
yyaxis left
h1 = plot(z_data, sphere_ratio, 'or', 'MarkerSize',4);
ylabel('\sigma_\theta / \sigma_\phi')
ylim([0 2.35])
grid on
hold on

% ---- RIGHT: geometry ----
yyaxis right
h2 = plot(z_geom_scaled, r_geom, 'k-', 'LineWidth',1.5);
hold on

% ---- Curvature ----
h3 = plot(z_geom_scaled, kappa_scaled, 'b--', 'LineWidth',1.5);

ylabel('Radius r (mm) / Scaled curvature')

xlabel('Longitudinal coordinate z (mm)')
title('Stress ratio variation with super-ellipsoid geometry and curvature')
xlim([0 930])

legend([h1 h2 h3], ...
       {'Stress ratio  \sigma_\theta / \sigma_\phi', ...
        'Super-ellipsoid profile', ...
        'Curvature (scaled)'}, ...
       'Location','northoutside', ...
       'Orientation','horizontal');

%% ---- Stress ratio from FE ----
Sphere_ratio = Sphere(:,4)./Sphere(:,3);
z_data = Sphere(:,2);

% ---- Super-ellipse parameters ----
a = 241.09;    % radial semi-axis (mm)
c = 241.09;    % longitudinal semi-axis (mm)
n = 1;      % exponent

theta = linspace(0, pi/2, 2000);

% ---- Geometry ----
r_geom = a .* (cos(theta).^n);
z_geom = c .* (sin(theta).^n);

% ---- First derivatives ----
dr_dth = -a * n .* (cos(theta).^(n-1)) .* sin(theta);
dz_dth =  c * n .* (sin(theta).^(n-1)) .* cos(theta);

% ---- Second derivatives ----
d2r_dth2 = -a*n*((n-1)*(cos(theta).^(n-2)).*(-sin(theta).^2) ...
                 + (cos(theta).^(n-1)).*cos(theta));

d2z_dth2 =  c*n*((n-1)*(sin(theta).^(n-2)).*(cos(theta).^2) ...
                 - (sin(theta).^(n-1)).*sin(theta));

% ---- Curvature ----
kappa = abs(dr_dth .* d2z_dth2 - dz_dth .* d2r_dth2) ./ ...
        ( (dr_dth.^2 + dz_dth.^2).^(3/2) );

% ---- Scale geometry to match FE axis ----
z_geom_scaled = z_geom / max(z_geom) * max(z_data);

% ---- Normalize curvature for visualization ----
kappa_scaled = kappa / max(kappa) * max(r_geom);

figure;

% ---- LEFT: stress ratio ----
yyaxis left
h1 = plot(z_data, Sphere_ratio, 'or', 'MarkerSize',4);
ylabel('\sigma_\theta / \sigma_\phi')
ylim([-1 3])
grid on
hold on

% ---- RIGHT: geometry ----
yyaxis right
h2 = plot(z_geom_scaled, r_geom, 'k-', 'LineWidth',1.5);
hold on

% ---- Curvature ----
h3 = plot(z_geom_scaled, kappa_scaled, 'b--', 'LineWidth',1.5);

ylabel('Radius r (mm) / Curvature')

xlabel('Longitudinal coordinate z (mm)')
title('Stress ratio variation with spherical geometry and curvature')
xlim([0 235])

legend([h1 h2 h3], ...
       {'Stress ratio  \sigma_\theta / \sigma_\phi', ...
        'Sphere profile', ...
        'Curvature (scaled)'}, ...
       'Location','northoutside', ...
       'Orientation','horizontal');
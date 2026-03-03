% ---- Super-ellipse parameters (EDIT) ----
a = 103.59;    % radial semi-axis (mm)
c = 932.3;    % longitudinal semi-axis (mm)
n = 0.65;      % exponent

theta = linspace(0, pi/2, 1000);

r_geom = a .* (cos(theta).^n);
z_geom = c .* (sin(theta).^n);

Sphere = readmatrix('hoopvmeri.xlsx', 'Range', 'B:E');
super_15 = readmatrix('hoopvmeri.xlsx', 'Range', 'G:J');
super_75 = readmatrix('hoopvmeri.xlsx', 'Range', 'L:O');

Sphere_ratio = Sphere(:,4)./Sphere(:,3);
figure;
plot(Sphere(:,2), Sphere_ratio, 'ok', 'DisplayName', 'Sphere Ratio');
xlim([1 235])
ylim([-1 3])
ylabel('\sigma_\theta / \sigma_\phi')
xlabel('Longitudinal coordinate z (mm)')
title('Stress ratio vs longitudinal coordinate for Sphere')
grid

super_15_ratio = super_15(:,4)./super_15(:,3);
z_data = super_15(:,2);
z_geom_scaled = z_geom / max(z_geom) * max(z_data);
figure;
plot(super_15(:,2), super_15_ratio, 'ob', 'DisplayName', 'Sphere Ratio');
xlim([1 930])
ylim([0 2.3])
ylabel('\sigma_\theta / \sigma_\phi')
grid

super_75_ratio = super_75(:,4) ./ super_75(:,3);
z_data = super_75(:,2);
z_geom_scaled = z_geom / max(z_geom) * max(z_data);

yyaxis right
plot(z_geom_scaled, r_geom, 'k-', 'LineWidth',1.5);
ylabel('Radius r (mm)')

xlabel('Longitudinal coordinate z (mm)')
title('Stress ratio vs longitudinal coordinate with super-ellipse profile')
xlim([min(z_data) max(z_data)])

% ---- Plot ----
figure

yyaxis left
h1 = plot(z_data, super_75_ratio, 'or', 'MarkerSize',4);
ylabel('\sigma_\theta / \sigma_\phi')
ylim([-5 3.2])
grid on
hold on

yyaxis right
h2 = plot(z_geom_scaled, r_geom, 'k-', 'LineWidth',1.5);
ylabel('Radius r (mm)')

xlabel('Longitudinal coordinate z (mm)')
title('Stress ratio vs longitudinal coordinate with super-ellipse profile')
xlim([min(z_data) max(z_data)])
legend([h1 h2], ...
       {'Stress ratio  \sigma_\theta / \sigma_m', ...
        'Super-ellipse profile  r(z)'}, ...
       'Location','northoutside', ...
       'Orientation','horizontal');
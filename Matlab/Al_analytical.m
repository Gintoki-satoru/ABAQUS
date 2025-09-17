% Ellipsoid parameters
a = 150;
b = 100;
t = 2.5;
P = 0.1;

phi = linspace(0, pi/2, 100);

x = a * cos(phi);
y = b * sin(phi);

R1 = ((b^2 - a^2)*x.^2 + a^4).^(3/2)/(a^4*b);
R2 = ((b^2 - a^2)*x.^2 + a^4).^(1/2)/b;

s_meri = P*R2/(2*t);
s_theta = (P/(2*t))*(((b^2 - a^2)*x.^2 + a^4).^(1/2)/b).*(2 - a^4./((b^2 - a^2)*x.^2 + a^4));

Fem_x = readmatrix('Circumfrential.xlsx', 'Range', 'B:C');
Fem_y = readmatrix('Circumfrential.xlsx', 'Range', 'E:F');
Fem_x_meri = readmatrix('Circumfrential.xlsx', 'Range', 'H:I');

% figure;
% plot(x, s_theta, 'b', 'LineWidth', 2,  'DisplayName','Analytical results'); hold on;
% plot(Fem_x(:,1), Fem_x(:,2), 'ro', 'MarkerSize', 4, 'DisplayName','FEM results');
% xlabel('x');
% ylabel('\sigma_\theta (Circumferential Stress)');
% % ylim([2.8 3.3])
% grid on;
% legend show
% title('Circumferential Stress \sigma_\theta vs x for Sphere');

% figure;
% plot(s_theta, y, 'b', 'LineWidth', 2,  'DisplayName','Analytical results'); hold on;
% plot(Fem_y(:,2), Fem_y(:,1), 'ro', 'MarkerSize', 4, 'DisplayName','FEM results');
% ylabel('y');
% xlabel('\sigma_\theta (Circumferential Stress)');
% grid on;
% legend show
% title('Circumferential Stress \sigma_\theta vs y');



data = readmatrix('element_stress_ip_coords.csv', 'Range', 'C:H');  % entire columns
data = data(2:end,:);  % skip the first row (header)

x_f = data(:,1);  % r-coordinate
y = data(:,2);  % z-coordinate
S11 = data(:,3); % σ_rr
S22 = data(:,4); % σ_zz
S12 = data(:,6); % σ_rz (shear)

n = length(x_f);

sigma_meridional = zeros(n,1);
sigma_normal = zeros(n,1);
tau_local = zeros(n,1);

for i = 1:n
    % --- 1. Calculate parametric angle phi for this point ---
    phi = atan2(y(i)/b, x_f(i)/a); % parametric angle along ellipse
 
    % --- 2. Compute tangent direction unit vector ---
    % Tangent vector = (-a*sin(phi), b*cos(phi))
    tx = -a*sin(phi);
    tz =  b*cos(phi);
    tnorm = sqrt(tx^2 + tz^2);
    tx = tx/tnorm; tz = tz/tnorm;

    % --- 3. Compute angle θ of tangent wrt r-axis ---
    theta = atan2(tz, tx);  % angle of tangent vector

    % --- 4. Stress tensor in global (r,z) ---
    sigma = [S11(i), S12(i);
             S12(i), S22(i)];

    % --- 5. Rotation matrix ---
    R = [cos(theta) sin(theta);
        -sin(theta) cos(theta)];

    % --- 6. Transform stress tensor ---
    sigma_prime = R * sigma * R'; 

    % --- 7. Store results ---
    sigma_meridional(i) = sigma_prime(1,1); % along tangent
    sigma_normal(i)     = sigma_prime(2,2); % along normal
    tau_local(i)        = sigma_prime(1,2); % shear in local frame
end

% Save results
output = [x_f, y, sigma_meridional, sigma_normal, tau_local];
writematrix(output, 'transformed_stresses.csv');

fprintf('Done! Results saved to transformed_stresses.csv\n');

figure;
plot(x, s_meri, 'b', 'LineWidth', 2,  'DisplayName','Analytical results'); hold on;
plot(x_f, sigma_meridional,'ro', 'MarkerSize', 4, 'DisplayName','FEM results');
xlabel('x');
ylabel('\sigma_\phi (Meridional Stress)');
% ylim([2 5])
grid on;
legend show
title('Meridional Stress \sigma_\phi vs x');
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

figure;
plot(x, s_theta, 'b', 'LineWidth', 2,  'DisplayName','Analytical results'); hold on;
plot(Fem_x(:,1), Fem_x(:,2), 'ro', 'MarkerSize', 4, 'DisplayName','FEM results');
xlabel('x');
ylabel('\sigma_\theta (Circumferential Stress)');
grid on;
legend show
title('Circumferential Stress \sigma_\theta vs x');

figure;
plot(s_theta, y, 'b', 'LineWidth', 2,  'DisplayName','Analytical results'); hold on;
plot(Fem_y(:,2), Fem_y(:,1), 'ro', 'MarkerSize', 4, 'DisplayName','FEM results');
ylabel('y');
xlabel('\sigma_\theta (Circumferential Stress)');
grid on;
legend show
title('Circumferential Stress \sigma_\theta vs y');


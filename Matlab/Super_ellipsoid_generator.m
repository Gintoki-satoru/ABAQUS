% Super-ellipsoid parameters
 a = 0.5761
 b = 0.5761
 c = 2.8806
n1 = 0.5;
n2 = 0.5;
thick = 0.00402;
epsilon1 = 0.5;
epsilon2 = 0.5;

signed_power = @(base, exp) sign(base) .* (abs(base).^exp);

% Parameter grid
[u, v] = meshgrid(linspace(-pi, pi, 80), linspace(-pi/2, pi/2, 40));

% Super-ellipsoid equations
x = a * signed_power(cos(v), epsilon2) .* signed_power(cos(u), epsilon1);
y = b * signed_power(cos(v), epsilon2) .* signed_power(sin(u), epsilon1);
z = c * signed_power(sin(v), epsilon2);

% Plot
figure;
surf(z, x, y, 'FaceColor', 'texturemap', 'EdgeColor', 'texturemap');
axis equal;
xlabel('Z'); ylabel('X'); zlabel('Y');
title('Super-Ellipsoid');
grid on;
camlight; lighting gouraud;

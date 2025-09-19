% Super-ellipsoid parameters
a = 3;  % semi-axis along x
b = 3;  % semi-axis along y
c = 2;  % semi-axis along z

epsilon1 = 1;  % controls roundness in x-y plane
epsilon2 = 1;  % controls roundness in vertical direction

signed_power = @(base, exp) sign(base) .* (abs(base).^exp);

% Parameter grid
[u, v] = meshgrid(linspace(-pi, pi, 80), linspace(-pi/2, pi/2, 40));

% Super-ellipsoid equations
x = a * signed_power(cos(v), epsilon2) .* signed_power(cos(u), epsilon1);
y = b * signed_power(cos(v), epsilon2) .* signed_power(sin(u), epsilon1);
z = c * signed_power(sin(v), epsilon2);

% Plot
figure;
surf(x, y, z, 'FaceColor', 'magenta', 'EdgeColor', 'none');
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Super-Ellipsoid');
grid on;
camlight; lighting gouraud;

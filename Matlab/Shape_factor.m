% Super-ellipsoid parameters
a = 100;    % semi-axis (m) or (mm, be consistent)
b = 100;
c = 500;
n1 = 4;
n2 = 4;

t = 10;  % thickness (same units as a,b,c)
k = 150;    % thermal conductivity (W/m·K or W/mm·K consistently!)

% Create grid
N = 200; % resolution (increase for accuracy)
phi = linspace(0, pi/2, N);
theta = linspace(0, pi/2, N);
[PHI, THETA] = meshgrid(phi, theta);

% Helper sign power function
spow = @(x,p) sign(x).*abs(x).^p;

% Parametric surface
X = a * spow(cos(PHI), 2/n1) .* spow(cos(THETA), 2/n2);
Y = b * spow(cos(PHI), 2/n1) .* spow(sin(THETA), 2/n2);
Z = c * spow(sin(PHI), 2/n1);

% Compute partial derivatives numerically
[dX_dphi, dX_dtheta] = gradient(X, phi, theta);
[dY_dphi, dY_dtheta] = gradient(Y, phi, theta);
[dZ_dphi, dZ_dtheta] = gradient(Z, phi, theta);

% r_phi and r_theta vectors
r_phi = cat(3, dX_dphi, dY_dphi, dZ_dphi);
r_theta = cat(3, dX_dtheta, dY_dtheta, dZ_dtheta);

% Cross product magnitude = local area density
crossprod = cross(r_phi, r_theta, 3);
dA = sqrt(sum(crossprod.^2, 3));

% Numerical integration
S = sum(dA(:)) * ( (phi(2)-phi(1)) * (theta(2)-theta(1)) ) / t;

% Thermal resistance
R = 1 / (k * S);

fprintf('Shape factor S = %.6g\n', S);
fprintf('Thermal resistance R = %.6g\n', R);

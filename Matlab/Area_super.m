% Super-ellipsoid parameters
a = 100;    % semi-axis (m) or (mm, be consistent)
b = 100;
c = 500;
n1 = 0.5;
n2 = 0.5;

% Parameter grid resolution (increase for accuracy)
N = 5250;
phi = linspace(0, pi/2, N);
theta = linspace(0, pi/2, N);
[PHI, THETA] = meshgrid(phi, theta);

% Helper: signed power function
spow = @(x,p) sign(x).*abs(x).^p;

% Parametric coordinates of surface
X = a * spow(cos(PHI), n1) .* spow(cos(THETA), n2);
Y = b * spow(cos(PHI), n1) .* spow(sin(THETA), n2);
Z = c * spow(sin(PHI), n1);

% Partial derivatives (numerical)
[dX_dphi, dX_dtheta] = gradient(X, phi, theta);
[dY_dphi, dY_dtheta] = gradient(Y, phi, theta);
[dZ_dphi, dZ_dtheta] = gradient(Z, phi, theta);

% Tangent vectors
r_phi   = cat(3, dX_dphi,   dY_dphi,   dZ_dphi);
r_theta = cat(3, dX_dtheta, dY_dtheta, dZ_dtheta);

% Local area element = |r_phi x r_theta|
cross_prod = cross(r_phi, r_theta, 3);
dA = sqrt(sum(cross_prod.^2, 3));

% Integrate over parameter domain
dphi = phi(2)-phi(1);
dtheta = theta(2)-theta(1);
A1 = sum(dA(:)) * dphi * dtheta;

fprintf('Surface Area = %.6f (square units)\n', A1);

%%
% Super-ellipsoid parameters
a = 110;    % semi-axis (m) or (mm, be consistent)
b = 110;
c = 510;
n1 = 0.5;
n2 = 0.5;

% Parameter grid resolution (increase for accuracy)
N = 5250;
phi = linspace(0, pi/2, N);
theta = linspace(0, pi/2, N);
[PHI, THETA] = meshgrid(phi, theta);

% Helper: signed power function
spow = @(x,p) sign(x).*abs(x).^p;

% Parametric coordinates of surface
X = a * spow(cos(PHI), n1) .* spow(cos(THETA), n2);
Y = b * spow(cos(PHI), n1) .* spow(sin(THETA), n2);
Z = c * spow(sin(PHI), n1);

% Partial derivatives (numerical)
[dX_dphi, dX_dtheta] = gradient(X, phi, theta);
[dY_dphi, dY_dtheta] = gradient(Y, phi, theta);
[dZ_dphi, dZ_dtheta] = gradient(Z, phi, theta);

% Tangent vectors
r_phi   = cat(3, dX_dphi,   dY_dphi,   dZ_dphi);
r_theta = cat(3, dX_dtheta, dY_dtheta, dZ_dtheta);

% Local area element = |r_phi x r_theta|
cross_prod = cross(r_phi, r_theta, 3);
dA = sqrt(sum(cross_prod.^2, 3));

% Integrate over parameter domain
dphi = phi(2)-phi(1);
dtheta = theta(2)-theta(1);
A2 = sum(dA(:)) * dphi * dtheta;

fprintf('Surface Area = %.6f (square units)\n', A2);

Alm = (A2 - A1)/(log(A2/A1))
Agm = sqrt(A2*A1)
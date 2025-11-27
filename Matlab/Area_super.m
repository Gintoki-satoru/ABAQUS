% Super-ellipsoid parameters
a = 161;
b = 161;
c = 725;
n1 = 1;
n2 = 1;

% Parameter grid resolution (increase for accuracy)
N = 5250;
phi = linspace(-pi/2, pi/2, N);
theta = linspace(-pi, pi, N);
[PHI, THETA] = meshgrid(phi, theta);

spow = @(x,p) sign(x).*abs(x).^p;

% Parametric coordinates of surface
X = a * spow(cos(PHI), n1) .* spow(cos(THETA), n2);
Y = b * spow(cos(PHI), n1) .* spow(sin(THETA), n2);
Z = c * spow(sin(PHI), n1);

% Partial derivatives
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
a = 100;
b = 100;
c = 500;
n1 = 0.5;
n2 = 0.5;
thick = 5;
V_inner = superellipsoid_volume(a, b, c, n1, n2);

V_outer = superellipsoid_volume(a+thick, b+thick, c+thick, n1, n2);

V_material = V_outer - V_inner;
fprintf('Superellipsoid Volume = %.3f mm^3\n', V_material);

% Parameter grid resolution
N = 5000;
phi = linspace(-pi/2, pi/2, N);
theta = linspace(-pi, pi, N);
[PHI, THETA] = meshgrid(phi, theta);

% Helper: signed power function
spow = @(x,p) sign(x).*abs(x).^p;

% Parametric coordinates of surface
X = a * spow(cos(PHI), n1) .* spow(cos(THETA), n2);
Y = b * spow(cos(PHI), n1) .* spow(sin(THETA), n2);
Z = c * spow(sin(PHI), n1);

% Partial derivatives (numerical)
[dX_dtheta, dX_dphi] = gradient(X, theta, phi);
[dY_dtheta, dY_dphi] = gradient(Y, theta, phi);
[dZ_dtheta, dZ_dphi] = gradient(Z, theta, phi);

% Tangent vectors
r_phi   = cat(3, dX_dphi,   dY_dphi,   dZ_dphi);
r_theta = cat(3, dX_dtheta, dY_dtheta, dZ_dtheta);

% Local area element = |r_phi x r_theta|
cross_prod = cross(r_phi, r_theta, 3);
dA = sqrt(sum(cross_prod.^2, 3));

% Integrate over parameter domain
dphi = phi(2)-phi(1);
dtheta = theta(2)-theta(1);
A = sum(dA(:)) * dphi * dtheta;

fprintf('Surface Area = %.3f mm^2\n', A);

S_inf = 3.51*sqrt(A);
S_0 = A/thick;
ls = 2 * max([a, b, c]);
expr = 1.26 - (2 - sqrt(A)/ls) / (9*sqrt(1 - 4.79*(V_inner)^(2/3)/A));
n = max(expr, 1.0);

S = (S_0^n + S_inf^n)^(1/n);

function V = superellipsoid_volume(a, b, c, n1, n2)
    e1 =2/n1;
    e2 =2/n2;
    V = 8 * a * b * c * (gamma(1 + 1/e1)^2 * gamma(1 + 1/e2)) / ...
        (gamma(1 + 2/e1) * gamma(1 + (1/e2 + 2/e1)));
end
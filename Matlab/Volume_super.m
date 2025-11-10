function V = superellipsoid_volume(a, b, c, n1, n2)
% SUPERELLSIPSOID_VOLUME  Computes the volume of a superellipsoid.
%
%   Volume formula:
%   V = 8 * a * b * c * (Gamma(1 + 1/n1)^2 * Gamma(1 + 1/n2)) ...
%           / (Gamma(1 + 2/n1) * Gamma(1 + (1/n2 + 2/n1)));
    e1 =2/n1;
    e2 =2/n2;
    V = 8 * a * b * c * (gamma(1 + 1/e1)^2 * gamma(1 + 1/e2)) / ...
        (gamma(1 + 2/e1) * gamma(1 + (1/e2 + 2/e1)));
end
a = 105;
b = 105;
c = 505;
n1 = 0.5;
n2 = 0.5;
thick = 5;

V_inner = superellipsoid_volume(a, b, c, n1, n2);

V_outer = superellipsoid_volume(a+thick, b+thick, c+thick, n1, n2);

V_material = V_outer - V_inner;
fprintf('Superellipsoid Volume = %.3f mm^3\n', V_material);


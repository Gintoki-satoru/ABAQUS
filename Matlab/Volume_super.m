function V = superellipsoid_volume(a, b, c, n1, n2)
% SUPERELLSIPSOID_VOLUME  Computes the volume of a superellipsoid.
%
%   Volume formula:
%   V = 8 * a * b * c * (Gamma(1 + 1/n1)^2 * Gamma(1 + 1/n2)) ...
%           / (Gamma(1 + 2/n1) * Gamma(1 + (1/n2 + 2/n1)));

    V = 8 * a * b * c * (gamma(1 + 1/n1)^2 * gamma(1 + 1/n2)) / ...
        (gamma(1 + 2/n1) * gamma(1 + (1/n2 + 2/n1)));
end
a = 100;
b = 100;
c = 500;
n1 = 4;
n2 = 4;
thick = 10;

V_inner = superellipsoid_volume(a, b, c, n1, n2);

V_outer = superellipsoid_volume(a+thick, b+thick, c+thick, n1, n2);

V_material = V_outer - V_inner;
fprintf('Superellipsoid Volume = %.3f mm^3\n', V_material);


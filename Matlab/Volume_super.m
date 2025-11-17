function V = superellipsoid_volume(a, b, c, n1, n2)
    e1 =2/n1;
    e2 =2/n2;
    V = 8 * a * b * c * (gamma(1 + 1/e1)^2 * gamma(1 + 1/e2)) / ...
        (gamma(1 + 2/e1) * gamma(1 + (1/e2 + 2/e1)));
end

 a = 0.141;
 b = 0.141;
 c = 0.705;
n1 = 1;
n2 = 1;
thick = 0.00189;

V_inner = superellipsoid_volume(a, b, c, n1, n2);

V_outer = superellipsoid_volume(a+thick, b+thick, c+thick, n1, n2);

V_material = V_outer - V_inner;
fprintf('Superellipsoid Volume = %.6f m^3\n', V_material);


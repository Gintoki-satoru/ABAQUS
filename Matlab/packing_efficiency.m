function V = superellipsoid_volume(a, b, c, n1, n2)
    e1 =2/n1;
    e2 =2/n2;
    V = 8 * a * b * c * (gamma(1 + 1/e1)^2 * gamma(1 + 1/e2)) / ...
        (gamma(1 + 2/e1) * gamma(1 + (1/e2 + 2/e1)));
end

a = 319.1;
b = 70.9;
c = 354.6;
n1 = 0.5;
n2 = 0.5;
thick = 5;

V_inner = superellipsoid_volume(a, b, c, n1, n2);

V_cube = 8*a*b*c;

eff = 5.870308322012401e+07/V_cube*100;
fprintf('Packing Efficiency = %.3f \n', eff);
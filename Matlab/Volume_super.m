function V = superellipsoid_volume(a1, a2, a3, n1, n2)
    eps1 = n1;
    eps2 = n2;

    B = @(p,q) beta(p,q);

    V = 2*a1*a2*a3 * eps1*eps2 * ...
        B(eps1/2 + 1, eps1) * ...
        B(eps2/2, eps2/2);
end

 a = 114.24;
 b = 114.24;
 c = 799.71;
n1 = 0.7;
n2 = 1;
thick = 1.576;

V_inner = superellipsoid_volume(a, b, c, n1, n2);

V_outer = superellipsoid_volume(a+thick, b+thick, c+thick, n1, n2);

V_material = V_outer - V_inner;
fprintf('Superellipsoid inner Volume = %.6f m^3\n', V_inner);
fprintf('Superellipsoid material Volume = %.6f m^3\n', V_material);


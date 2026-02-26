function V = superellipsoid_volume(a1, a2, a3, n1, n2)
    eps1 = n1;   % controls z-related rounding
    eps2 = n2;   % controls xy-plane rounding

    % Beta function (MATLAB has beta(); otherwise use gamma relation)
    B = @(p,q) beta(p,q);
    % If you prefer no beta(), uncomment and use:
    % B = @(p,q) gamma(p).*gamma(q)./gamma(p+q);

    V = 2*a1*a2*a3 * eps1*eps2 * ...
        B(eps1/2 + 1, eps1) * ...
        B(eps2/2, eps2/2);
end

 a = 87.71;
 b = 87.71;
 c = 1227.89;
n1 = 0.4;
n2 = 1;
thick = 1.252;

V_inner = superellipsoid_volume(a, b, c, n1, n2);

V_outer = superellipsoid_volume(a+thick, b+thick, c+thick, n1, n2);

V_material = V_outer - V_inner;
fprintf('Superellipsoid Volume = %.6f m^3\n', V_material);


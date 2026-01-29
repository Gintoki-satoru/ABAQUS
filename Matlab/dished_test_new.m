a   = 114.24;
b   = 799.71;
n1   = 2/0.7;
h  = 2;
p0 = 1.0;      % MPa 
N  = 2000;

out = head_superellipse_stress(a, b, h, p0, n1, N);

% filter out NaN/Inf at rho=0 or rho=1
ok = isfinite(out.sigma_theta) & isfinite(out.sigma_phi) & isfinite(out.sigma_VM);

r  = out.r(ok);
st = out.sigma_theta(ok);
sp = out.sigma_phi(ok);
sv = out.sigma_VM(ok);

figure;
plot(r, st, 'LineWidth', 1.5); grid on;
xlabel('r (mm)'); ylabel('\sigma_\theta (MPa)');
title('Meridional membrane stress \sigma_\theta');

figure;
plot(r, sp, 'LineWidth', 1.5); grid on;
xlabel('r (mm)'); ylabel('\sigma_\phi (MPa)');
title('Hoop membrane stress \sigma_\phi');

figure;
plot(r, sv, 'LineWidth', 1.5); grid on;
xlabel('r (mm)'); ylabel('\sigma_{VM} (MPa)');
title('Von Mises stress (membrane, plane stress)');



function out = head_superellipse_stress(a, b, h, p0, n1, N)

    if nargin < 6, N = 200; end

    beta = b / a;
    rho  = linspace(0.05, 0.99, N).';
    r    = rho * a;

    % Dimensionless meridian: x/a = g(rho) = beta * [1 - (1 - rho^n1)^(1/n1)]
    rho_n = rho.^n1;
    one_m = 1 - rho_n;

    g  = beta .* ( 1 - one_m.^(1./n1) );

    % Derivatives g'(rho), g''(rho) for general Lamé exponent n1
    % From earlier derivation:
    % g'(rho)  = beta * rho^(n1-1) * (1 - rho^n1)^(1/n1 - 1)
    % g''(rho) = beta * (n1-1) * rho^(n1-2) * (1 - rho^n1)^(1/n1 - 2)
    g1 = beta .* rho.^(n1-1) .* one_m.^(1./n1 - 1);
    g2 = beta .* (n1-1) .* rho.^(n1-2) .* one_m.^(1./n1 - 2);

    % Curvatures
    A  = 1 + g1.^2;
    R1_over_a = (A).^(3/2) ./ g2;
    R2_over_a = rho .* sqrt(A) ./ g1;

    R1 = a * R1_over_a;
    R2 = a * R2_over_a;

    % Membrane stresses
    sigma_theta = p0 * R2 ./ (2*h);
    sigma_phi   = (2 - R2./R1) .* sigma_theta;
    sigma_VM    = sqrt( sigma_theta.^2 - sigma_theta.*sigma_phi + sigma_phi.^2 );

    % Physical x coordinate
    x = g * a;

    out.r         = r;
    out.x         = x;
    out.rho       = rho;
    out.R1        = R1;
    out.R2        = R2;
    out.sigma_theta = sigma_theta;
    out.sigma_phi   = sigma_phi;
    out.sigma_VM    = sigma_VM;
end

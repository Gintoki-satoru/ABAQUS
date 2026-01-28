% ---- Inputs ----
a   = 241.09;   % radius at cylinder junction (mm)
b   = 241.09;   % head depth (mm)
n   = 2.5;        % super-ellipse exponent (same power in both terms)
h   = 2;    % wall thickness (mm)
p0  = 1.0;      % pressure (MPa if you want stresses in MPa) (consistent units!)

N   = 2000;     % number of points along meridian

% ---- Dimensionless coordinate ----
beta = b/a;
rho  = linspace(1e-6, 1-1e-6, N);   % avoid rho=0 and rho=1 singularities
r    = a*rho;

% ---- Geometry: slope ----
% A(rho) = tan(theta) = dx/dr
A = beta .* (rho.^(n-1)) .* (1 - rho.^n).^((1/n) - 1);

% Helpful factors
onePlusA2 = 1 + A.^2;
sqrt1pA2  = sqrt(onePlusA2);

% ---- Curvatures / radii ----
% R2/a = rho * sqrt(1+A^2) / A
R2 = a .* (rho .* sqrt1pA2 ./ A);

% R1/a = (1+A^2)^(3/2) / [ beta*(n-1)*rho^(n-2)*(1-rho^n)^(1/n - 2) ]
den = beta*(n-1) .* (rho.^(n-2)) .* (1 - rho.^n).^((1/n) - 2);
R1 = a .* ((onePlusA2).^(3/2) ./ den);

% Ratio
R2_over_R1 = R2 ./ R1;
% Closed form:
R2_over_R1_closed = (n-1) ./ ((1 - rho.^n) .* (1 + A.^2));

% ---- Membrane stresses ----
sigma_theta = (R2 ./ (2*h)) .* p0;                 
sigma_phi   = (2 - R2_over_R1) .* sigma_theta;     % hoop

% ---- Optional: x(z) meridian coordinate (for plotting shape) ----
% x/a = beta * (1 - (1 - rho^n)^(1/n))
x = a * beta .* (1 - (1 - rho.^n).^(1/n));  % x measured from apex (x=0 at rho=0)

% ---- Quick checks ----
max_rel_err = max(abs(R2_over_R1 - R2_over_R1_closed) ./ max(1e-12, abs(R2_over_R1_closed)));
fprintf('Max relative error (R2/R1 numeric vs closed-form): %.3e\n', max_rel_err);

% ---- Plots ----
figure; plot(r, R1, 'LineWidth', 1.5); grid on;
xlabel('r (mm)'); ylabel('R_1 (mm)'); title('Meridional radius of curvature R_1');

figure; plot(r, R2, 'LineWidth', 1.5); grid on;
xlabel('r (mm)'); ylabel('R_2 (mm)'); title('Hoop radius of curvature R_2');

figure; plot(r, sigma_theta, 'LineWidth', 1.5); grid on;
xlabel('r (mm)'); ylabel('\sigma_\theta (MPa)'); title('Meridional membrane stress');

figure; plot(r, sigma_phi, 'LineWidth', 1.5); grid on;
xlabel('r (mm)'); ylabel('\sigma_\phi (MPa)'); title('Hoop membrane stress');

% Shape plot (x vs r)
figure; plot(x, r, 'LineWidth', 1.5); grid on; axis equal;
xlabel('x (mm)'); ylabel('r (mm)'); title('Meridian profile (apex to cylinder junction)');

%% ---- Outputs in a struct (handy for saving) ----
out.r = r; out.rho = rho; out.x = x;
out.R1 = R1; out.R2 = R2;
out.sigma_theta = sigma_theta; out.sigma_phi = sigma_phi;
assignin('base','out',out);

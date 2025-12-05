a = 124.34;     % mm
beta = 0.7;
b = 621.7;
k = 0.3;
h   = 1.67;   % [mm]
p0  = 1;     % [Pa] = 1 MPa
N   = 200;

[rho, s_th, s_ph, s_VM] = ellipsoidal_head_vonmises(a,b,h,p0,k,N);

% quick plot of von Mises vs parallel radius r = rho*a
figure
r = rho * a;
plot(r, s_VM); xlabel('r [m]'); ylabel('\sigma_{VM} [MPa]');
grid on;

plot_head_shape(a, b, k);



function [rho, sigma_theta, sigma_phi, sigma_VM] = ...
    ellipsoidal_head_vonmises(a, b, h, p0, k, N)
% ELLIPSOIDAL_HEAD_VONMISES
% Analytical membrane Huber–Mises (von Mises) stress in a
% generalised ellipsoidal dished head, after Magnucki & Jasion (2022).
%
% Geometry / loading:
%   a  - cylinder/head radius   [length]
%   b  - head depth             [length]
%   h  - shell thickness        [length]
%   p0 - pressure magnitude     [stress] (internal or external, sign irrelevant for VM)
%   k  - shape parameter 0 < k <= 0.5   (Eq. (2) in the paper)
%   N  - number of points along meridian (optional, default = 200)
%
% Outputs (vectors of length N):
%   rho         - dimensionless parallel radius r/a, in [0,1]
%   sigma_theta - meridional (longitudinal) membrane stress [same units as p0]
%   sigma_phi   - circumferential (hoop) membrane stress   [same units as p0]
%   sigma_VM    - Huber–Mises equivalent stress            [same units as p0]
%
% Reference:
% K. Magnucki, P. Jasion, "Strength of a cylindrical pressure vessel
% with individual ellipsoidal dished heads", Int. J. Pressure Vessels
% and Piping, 199 (2022) 104751.

    if nargin < 6
        N = 200;
    end

    % Dimensionless parameters
    beta = b / a;                % beta = b/a
    rho  = linspace(0, 1, N).';  % 0 <= rho = r/a <= 1

    % Common expressions from Eqs. (12)–(13)
    % A = (2 k beta)^2 * rho^2 + (1 - rho^2)^(2(1-k))
    two_k_beta = 2 * k * beta;
    rho2       = rho.^2;
    one_minus_rho2 = 1 - rho2;

    A = (two_k_beta^2) * rho2 + one_minus_rho2.^(2*(1 - k));

    % ---- Meridional (longitudinal) membrane stress sigma_theta (Eq. (12)) ----
    % sigma_theta(rho) = [1/(4 k beta)] * sqrt(A) * (a/h) * p0
    sigma_theta = (1 ./ (4 * k * beta)) .* sqrt(A) .* (a / h) * p0;

    % ---- Circumferential (hoop) membrane stress sigma_phi (Eq. (13)) ----
    % factor inside brackets:
    %   2 - { [1 + (1 - 2k) rho^2] * (1 - rho^2)^(1 - 2k) / A }
    num_f = (1 + (1 - 2*k) * rho2) .* one_minus_rho2.^(1 - 2*k);
    fac   = 2 - num_f ./ A;

    sigma_phi = fac .* sigma_theta;

    % ---- Huber–Mises equivalent stress (Eq. (14)) ----
    % sigma_VM = sqrt( sigma_theta^2 - sigma_theta*sigma_phi + sigma_phi^2 )
    sigma_VM = sqrt( sigma_theta.^2 - sigma_theta .* sigma_phi + sigma_phi.^2 );
end

function plot_head_shape(a, b, k, N)
% PLOT_HEAD_SHAPE  Plot meridian of generalized ellipsoidal dished head
% Based on Eq. (2) from Magnucki & Jasion (2022).
%
% a = head radius at junction [mm or m]
% b = head depth             [same units]
% k = shape parameter (0 < k <= 0.5)
% N = number of sample points (optional)

    if nargin < 4
        N = 300;
    end

    % r goes from tip (r=0) to cylinder junction (r=a)
    r = linspace(0, a, N).';

    % Eq. (2): x / b = 1 - [1 - (r/a)^2]^k
    x = b * (1 - (1 - (r/a).^2).^k);

    % Plot
    figure; hold on; grid on;
    plot(x, r, 'LineWidth', 2);
    % set(gca, 'YDir','reverse');   % optional: match paper orientation
    xlabel('x'); ylabel('r');
    title(sprintf('Generalized Ellipsoidal Head Shape (a=%.1f, b=%.1f, k=%.3f)', a, b, k));
    axis equal;

end

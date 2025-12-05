%% MAIN SCRIPT: Compare paper head vs. equivalent superellipse head

clear; clc; close all;

% Geometry & loading (consistent with the paper: mm & MPa)
a   = 126.01;      % [mm] cylinder/head radius
b   = 630.07; % [mm] head depth
h   = 1.693;        % [mm] thickness
p0  = 1;        % [MPa] internal pressure
k   = 0.65/2;      
N   = 400;      % points along meridian

% 1) Fit n1 so that superellipse meridian best matches paper meridian
n1_fit = fit_superellipse_n1(a, b, k);
fprintf('Fitted superellipse exponent n1 = %.4f\n', n1_fit);

% 2) Compute stresses for both geometries
out_paper = head_paper_stress(a, b, h, p0, k, N);
out_super = head_superellipse_stress(a, b, h, p0, 2/0.65, N);

r = out_paper.r;   % same sampling r for both

a3D = a;
c3D = b;
n1 = n1_fit;


function [r_barr, x_barr] = barr_meridian_axisymmetric(a3D, c3D, n1, N)

    if nargin < 4
        N = 400;
    end

    % signed power
    spow = @(u,p) sign(u).*abs(u).^p;

    % φ sweeps the meridian
    phi = linspace(-pi/2, pi/2, N).';

    % AXISYMMETRIC BARR MERIDIAN (n2 = 1)
    x_barr = a3D .* spow(cos(phi), 2/n1);
    r_barr = c3D .* spow(sin(phi), 2/n1);

    % only positive radius for the r-x plot
    r_barr = abs(r_barr);
end


%% 3) Plot shape comparison
% Compute Barr meridian cross-section
[r_barr, x_barr] = barr_meridian_axisymmetric(a3D, c3D, n1_fit, 400);


% Plot all curves together
figure; hold on; grid on; axis equal;

plot(out_paper.r,   out_paper.x, 'b-',  'LineWidth', 1.8);
plot(out_super.r,   out_super.x, 'r--', 'LineWidth', 1.8);
plot(x_barr, r_barr, 'k:', 'LineWidth', 1.8);

xlabel('r [mm]');
ylabel('x [mm]');

title(sprintf('Meridian comparison: paper k = %.3f, superellipse n1 = %.3f', ...
               k, n1_fit));

legend('Paper head', 'Equivalent superellipse (Lamé)', ...
       '3D Barr superellipsoid meridian', ...
       'Location', 'Best');


%% 4) Plot von Mises stress comparison
figure; hold on; grid on;
plot(r, out_paper.sigma_VM, 'b-', 'LineWidth', 1.8);
plot(r, out_super.sigma_VM, 'r--', 'LineWidth', 1.8);
xlabel('r [mm]');
ylabel('\sigma_{VM} [MPa]');
title('Von Mises stress along meridian');
legend('Paper head','Equivalent superellipse','Location','Best');


%% ========================================================================
%  Fit n1: best superellipse (Lamé) meridian to match paper head meridian
% ========================================================================
function n1_fit = fit_superellipse_n1(a, b, k)

    % Dimensionless grid
    Nfit = 300;
    rho  = linspace(0, 1, Nfit).';
    
    % Paper meridian (dimensionless X = x/b)
    X_paper = 1 - (1 - rho.^2).^k;

    % Objective: least-squares error between X_se(rho; n1) and X_paper
    fun = @(n1) sum( ( X_superellipse(rho, n1) - X_paper ).^2 );
    
    n1_init = 2.0;  % start near ellipse/circle
    n1_fit  = fminsearch(fun, n1_init, optimset('Display','off'));
    n1_fit  = abs(n1_fit);  % enforce positivity
end

% Superellipse (Lamé) meridian in non-dimensional form:
% |X - 1|^n1 + |rho|^n1 = 1  =>  X = 1 - (1 - rho^n1)^(1/n1)
function X = X_superellipse(rho, n1)
    X = 1 - (1 - rho.^n1).^(1./n1);
end


%% ========================================================================
%  Paper head stresses: x/b = 1 - [1 - (r/a)^2]^k  (Magnucki & Jasion)
% ========================================================================
function out = head_paper_stress(a, b, h, p0, k, N)

    if nargin < 6, N = 200; end

    beta = b / a;
    rho  = linspace(0, 1, N).';       % rho = r/a
    r    = rho * a;                   % [mm]

    % Dimensionless meridian: g(rho) = x/a
    g  = beta * ( 1 - (1 - rho.^2).^k );
    
    % Derivatives g'(rho), g''(rho) analytically
    % g'(rho)  = 2*k*beta*rho*(1 - rho^2)^(k-1)
    % g''(rho) = 2*k*beta*(1 - rho^2)^(k-2) * [1 - (2k-1)*rho^2]
    one_m = 1 - rho.^2;
    g1 = 2 * k * beta .* rho .* one_m.^(k-1);
    g2 = 2 * k * beta .* one_m.^(k-2) .* (1 - (2*k-1).*rho.^2);

    % Curvatures (using x/a = g, x' = g1, x'' = g2/a):
    % R1/a = (1 + g1^2)^(3/2) / g2
    A  = 1 + g1.^2;
    R1_over_a = (A).^(3/2) ./ g2;
    R2_over_a = rho .* sqrt(A) ./ g1;

    R1 = a * R1_over_a;
    R2 = a * R2_over_a;

    % Membrane stresses
    sigma_theta = p0 * R2 ./ (2*h);          % [MPa]
    sigma_phi   = (2 - R2./R1) .* sigma_theta;
    sigma_VM    = sqrt( sigma_theta.^2 - sigma_theta.*sigma_phi + sigma_phi.^2 );

    % Physical x coordinate
    x = g * a;   % since g = x/a

    out.r         = r;
    out.x         = x;
    out.rho       = rho;
    out.R1        = R1;
    out.R2        = R2;
    out.sigma_theta = sigma_theta;
    out.sigma_phi   = sigma_phi;
    out.sigma_VM    = sigma_VM;
end


%% ========================================================================
%  Superellipse (Lamé) head stresses: x/b = 1 - [1 - (r/a)^n1]^(1/n1)
% ========================================================================
function out = head_superellipse_stress(a, b, h, p0, n1, N)

    if nargin < 6, N = 200; end

    beta = b / a;
    rho  = linspace(0, 1, N).';
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

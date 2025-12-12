% USER INPUTS
a   = 241.09;
b   = 241.09;
c   = 241.09;
t   = 2.372;
n1  = 1;
n2  = 1;

% CONSTANT PARAMETERS
t_ins   = 16;
t_outer = 2;

k_liner = 0.0306;        % W/mm-K
k_outer = 0.123;        % W/mm-K
k_ins   = 3.03e-08;      % W/mm-K

T_air = 300;
T_LH2 = 20;

hc = 10 / 1e6;           % W/(mm^2-K)

% FUNCTION HANDLES FOR SUPERELLIPSOID AREA & VOLUME
superellipsoid_area = @(a,b,c,n1,n2) area_superellipsoid(a,b,c,n1,n2);
superellipsoid_volume = @(a,b,c,n1,n2) volume_superellipsoid(a,b,c,n1,n2);

% 1. SURFACE AREAS
A_liner = superellipsoid_area(a, b, c, n1, n2);
A_ins   = superellipsoid_area(a+t, b+t, c+t, n1, n2);
A_outer = superellipsoid_area(a+t+t_ins, b+t+t_ins, c+t+t_ins, n1, n2);
A_outer_total = superellipsoid_area(a+t+t_ins+t_outer, b+t+t_ins+t_outer, c+t+t_ins+t_outer, n1, n2);

% 2. VOLUMES
V_inner = superellipsoid_volume(a+t, b+t, c+t, n1, n2) ...
        - superellipsoid_volume(a, b, c, n1, n2);

V_ins = superellipsoid_volume(a+t+t_ins, b+t+t_ins, c+t+t_ins, n1, n2) ...
      - superellipsoid_volume(a+t, b+t, c+t, n1, n2);

V_outer = superellipsoid_volume(a+t+t_ins+t_outer, b+t+t_ins+t_outer, c+t+t_ins+t_outer, n1, n2) ...
        - superellipsoid_volume(a+t+t_ins, b+t+t_ins, c+t+t_ins, n1, n2);

% 3. SHAPE FACTORS
S_liner = shape_factor(A_liner, V_inner, t, a, b, c);
S_ins   = shape_factor(A_ins,   V_ins,  t_ins, a, b, c);
S_outer = shape_factor(A_outer, V_outer, t_outer, a, b, c);

% 4. HEAT TRANSFER SOLUTION
[Q_total, T3, T2, T1, h_eq] = equivalent_heat_coeff(...
    T_air, T_LH2, A_outer_total, hc,...
    S_ins, k_ins, S_liner, k_liner,...
    S_outer, k_outer, A_ins);
Bor = Q_total*8640000/(445600*4.1);

% OUTPUT RESULTS
fprintf('\n===================== RESULTS =====================\n');
fprintf('Input Geometry: a=%.2f, b=%.2f, c=%.2f, t=%.2f, n1=%.2f, n2=%.2f\n',a,b,c,t,n1,n2);
fprintf('-----------------------------------------------------\n');
fprintf('Total Heat Load Q_total = %.6f W\n', Q_total);
fprintf('BOR = %.6f\n', Bor);
fprintf('=====================================================\n');



% SUPERELLIPSOID AREA FUNCTION
function A = area_superellipsoid(a,b,c,n1,n2)
    N = 500;
    phi   = linspace(-pi/2, pi/2, N);
    theta = linspace(-pi, pi, N);
    [PHI, THETA] = meshgrid(phi, theta);
    sp = @(x,p) sign(x).*abs(x).^p;
    X = a*sp(cos(PHI),n1).*sp(cos(THETA),n2);
    Y = b*sp(cos(PHI),n1).*sp(sin(THETA),n2);
    Z = c*sp(sin(PHI),n1);
    [dX_dphi, dX_dtheta] = gradient(X, phi, theta);
    [dY_dphi, dY_dtheta] = gradient(Y, phi, theta);
    [dZ_dphi, dZ_dtheta] = gradient(Z, phi, theta);
    rphi   = cat(3, dX_dphi, dY_dphi, dZ_dphi);
    rtheta = cat(3, dX_dtheta, dY_dtheta, dZ_dtheta);
    cp = cross(rphi, rtheta, 3);
    dA = sqrt(sum(cp.^2,3));
    A = sum(dA,'all')*(phi(2)-phi(1))*(theta(2)-theta(1));
end


% SUPERELLIPSOID VOLUME FUNCTION
function V = volume_superellipsoid(a,b,c,n1,n2)
    e1 = 2/n1;
    e2 = 2/n2;
    V = 8*a*b*c*(gamma(1+1/e1))^2*gamma(1+1/e2) ...
        /(gamma(1+2/e1)*gamma(1+(1/e2+2/e1)));
end


% SHAPE FACTOR
function S = shape_factor(A,Vinner,t,a,b,c)
    S_inf = 3.51*sqrt(A);
    S0 = A/t;
    ls = 2*max([a,b,c]);
    expr = 1.26-(2-sqrt(A)/ls)/(9*sqrt(1-4.79*(Vinner^(2/3))/A));
    n = max(expr,1);
    S = (S0^n + S_inf^n)^(1/n);
end


% EQUIVALENT HEAT COEFFICIENT
function [Q_total,T3,T2,T1,h_eq] = equivalent_heat_coeff(T_air,T_LH2,A_outer,hc,S_ins,k_ins,S_liner,k_liner,S_outer,k_outer,A_outer_liner)
    R2 = 1/(A_outer*hc);
    R1 = 1/(k_ins*S_ins);
    Rliner = 1/(k_liner*S_liner);
    Router = 1/(k_outer*S_outer);
    Aeq = [R2 1 0 0;
           Router -1 1 0;
           R1 0 -1 1;
           Rliner 0 0 -1];
    beq = [T_air; 0; 0; -T_LH2];
    x = Aeq\beq;
    Q = x(1);
    T3 = x(2); T2 = x(3); T1 = x(4);
    Q_total = 1.13 * Q * 1.3;
    h_eq = Q_total / ((T_air - T_LH2)*A_outer_liner);
end

%%
a   = 241.09;
b   = 241.09;
c   = 241.09;
t   = 2.372;
n1  = 1;
n2  = 1;

% CONSTANT PARAMETERS
t_ins   = 16;
t_outer = 2;

k_liner = 0.0306;        % W/mm-K
k_outer = 0.123;        % W/mm-K
k_ins   = 3.03e-08;      % W/mm-K

T_air = 300;
T_LH2 = 20;

hc = 10 / 1e6;
P = 1e-2;              % Pa (given)

% 2) k_e(P) for each insulation material (from figure)
%    Units in the paper: mW/(m·K)
%    YOU must fill in the exact formulas from the paper.

% Example structure (replace "..." with real expressions):
% ke_fabric_mW = 0.13 + (0.1*P)./(P + 0.02) + (11.9*P)./(P + 23.5);
% etc.

% --- FILL THESE FROM THE PAPER ---
ke_fabric_mW  = 0.13+ (0.1*P)/(P+0.02)+(11.9*P)/(P+23.5);  % 15-layer fabric/foil
ke_MLI40_mW   = 0.07+ (1.9*P)/(P+8.8)+(13.4*P)/(P+105);  % 40-layer MLI
ke_foil_mW    = 0.06+ (14.8*P)/(P+179.3)+(9.4*P)/(P+120165);  % MLI foil paper
ke_FG_mW      = 1.95+ (14.5*P)/(P+22.2)+(11.8*P)/(P+23499);  % fibre glass
ke_GB_mW      = 0.68+ (1.3*P)/(P+39.8)+(23.7*P)/(P+398.4);  % glass bubbles
ke_MLInet_mW  = 0.03+ (9.8*P)/(P+28.4)+(8.7*P)/(P+16528);  % MLI Mylar-net
ke_AB_mW      = 1.7+ (3.6*P)/(P+33.8)+(10.3*P)/(P+70535);  % aerogel blanket
ke_SOF_mW     = 7.3+ (12.6*P)/(P+1);  % spray-on foam

% Put them in a vector
ke_mW = [ke_fabric_mW, ke_MLI40_mW, ke_foil_mW, ...
         ke_FG_mW, ke_GB_mW, ke_MLInet_mW, ke_AB_mW, ke_SOF_mW];

matNames = { ...
    'Fabric/Foil', ...
    '40-layer MLI', ...
    'MLI Foil Paper', ...
    'Fiber Glass', ...
    'Glass Bubbles', ...
    'MLI Mylar-Net', ...
    'Aerogel Blanket', ...
    'Spray-on Foam' };

% ---------------------------------------------------------
% Convert ke from mW/(m·K) to W/(mm·K)
% 1 mW/m·K = 1e-3 W/m·K = 1e-3 / 1000 W/mm·K = 1e-6 W/mm·K
% ---------------------------------------------------------
ke_W_per_mmK = ke_mW * 1e-6;  % this will be used as k_ins per material

% 3) Geometry-based quantities that do NOT depend on material
% Assuming you already have these function on your path
%   area_superellipsoid(a,b,c,n1,n2)
%   volume_superellipsoid(a,b,c,n1,n2)
%   shape_factor(A,V,t,a,b,c)
%   equivalent_heat_coeff(...)

A_liner = area_superellipsoid(a, b, c, n1, n2);
A_ins   = area_superellipsoid(a+t, b+t, c+t, n1, n2);
A_outer = area_superellipsoid(a+t+t_ins, b+t+t_ins, c+t+t_ins, n1, n2);
A_outer_total = area_superellipsoid(a+t+t_ins+t_outer, ...
                                    b+t+t_ins+t_outer, ...
                                    c+t+t_ins+t_outer, n1, n2);

V_inner = volume_superellipsoid(a+t, b+t, c+t, n1, n2) ...
        - volume_superellipsoid(a, b, c, n1, n2);

V_ins = volume_superellipsoid(a+t+t_ins, b+t+t_ins, c+t+t_ins, n1, n2) ...
      - volume_superellipsoid(a+t, b+t, c+t, n1, n2);

V_outer = volume_superellipsoid(a+t+t_ins+t_outer, ...
                                b+t+t_ins+t_outer, ...
                                c+t+t_ins+t_outer, n1, n2) ...
        - volume_superellipsoid(a+t+t_ins, ...
                                b+t+t_ins, ...
                                c+t+t_ins, n1, n2);

S_liner = shape_factor(A_liner, V_inner, t,       a, b, c);
S_ins   = shape_factor(A_ins,   V_ins,  t_ins,   a, b, c);
S_outer = shape_factor(A_outer, V_outer, t_outer, a, b, c);

% 4) Loop over materials, compute Q_total and BOR
nMat = numel(ke_W_per_mmK);
Q_total = zeros(1, nMat);
Bor     = zeros(1, nMat);

for i = 1:nMat
    k_ins_i = ke_W_per_mmK(i);  % effective insulation conductivity

    [Q_i, T3_i, T2_i, T1_i, h_eq_i] = equivalent_heat_coeff( ...
        T_air, T_LH2, A_outer_total, hc, ...
        S_ins, k_ins_i, ...
        S_liner, k_liner, ...
        S_outer, k_outer, ...
        A_ins);

    Q_total(i) = Q_i;

    % BOR = Q_total*8640000/(445600*4.1);
    Bor(i) = Q_i * 8640000 / (445600 * 4.1);
end

% 5) Plot BOR for each material
figure;
b = bar(Bor);
grid on;
set(gca,'XTick',1:nMat,'XTickLabel',matNames,'XTickLabelRotation',35);
ylabel('BOR % per day');
title(sprintf('Boil off rate(BOR) per day for different insulation materials at P = %.1e Pa', P));

% --- Add value labels above each bar ---
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = string(round(b.YData, 3));  % round to 3 decimals (adjust if needed)

text(xtips, ytips, labels, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontSize',10);

%% BOR vs slenderness

function plot_BOR_vs_slenderness(csvfile)

    % Load data
    T = readtable(csvfile);

    % --- Restrict to n-range 0.6 to 1.0 ---
    T = T(T.n1 >= 0.6 & T.n1 <= 1.0 & ...
          T.n2 >= 0.6 & T.n2 <= 1.0, :);

    % Extract
    a = T.a_mm;
    b = T.b_mm;
    c = T.c_mm;

    % Only keep axisymmetric shapes (a = b)
    mask_ab = abs(a - b) < 1e-6;
    T = T(mask_ab, :);

    % Compute slenderness λ = c/a
    lambda = T.c_mm ./ T.a_mm;

    % Remove undesired λ < 1 values
    T = T(lambda >= 1, :);
    lambda = lambda(lambda >= 1);

    % ----------------------------------------------------
    % PLOT: n2 = 1, n1 varies, λ = [1,3,5,7,9]
    % ----------------------------------------------------

    slender_target_1 = [1 3 5 7 9];
    n2_fixed = 1;
    n1_values = unique(T.n1);

    figure; hold on;

    for n1_val = n1_values'

        mask = abs(T.n1 - n1_val) < 1e-6 & abs(T.n2 - n2_fixed) < 1e-6;

        if ~any(mask), continue; end

        lam = lambda(mask);
        BOR = T.BOR(mask);   % <------ CHANGE HERE

        % Select exact slenderness targets (rounded)
        keep = ismember(round(lam), slender_target_1);
        lam  = lam(keep);
        BOR  = BOR(keep);

        if numel(lam) < 2, continue; end

        % Sort before interpolation
        [lam_sorted, idx] = sort(lam);
        BOR_sorted = BOR(idx);

        % Spline interpolation for smooth curve
        lam_fine = linspace(min(lam_sorted), max(lam_sorted), 200);
        BOR_fine = spline(lam_sorted, BOR_sorted, lam_fine);

        plot(lam_fine, BOR_fine, 'LineWidth', 1.5);
        plot(lam_sorted, BOR_sorted, 'o', 'DisplayName', sprintf('n1 = %.2f', n1_val));
    end

    xlabel('Slenderness ratio \lambda = c/a')
    ylabel('Boil-off rate BOR (%/day)')
    title('BOR vs Slenderness ratio (n2 = 1, n1 decreasing)')
    legend show; grid on;
    hold off;
end

plot_BOR_vs_slenderness("cleaned_file.csv")
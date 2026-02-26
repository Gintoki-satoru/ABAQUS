% Known temperatures
T_air = 300;     % K
T_LH2 = 20;      % K

% OUTER CONVECTION
A  = 1176179.901067;       % mm^2 (outer surface)
hc = 10;             % W/m^2K
hc_mm = hc/1e6;      % W/mm^2-K
R2 = 1 / (A * hc_mm);

% INSULATION
k_eff_mW = 0.0303;        % mW/mK
k_eff = k_eff_mW * 1e-6;  % W/mmK
S_ins = 67003.022;  % mm
G1 = k_eff * S_ins;
R1 = 1 / G1;

% LINER
k_al = 0.0306;        % W/mmK
S_liner = 502540.155;  % mm
G_liner = k_al * S_liner;
Rliner = 1 / G_liner;

% OUTER WALL METAL RESISTANCE
k_outer = 0.0306;       % W/mmK
S_outer = 582749.652;       % mm
Gouter  = k_outer * S_outer;
Router  = 1 / Gouter;

% Unknowns vector: x = [Q; T3; T2; T1]
% Equations:
% 1) R2*Q + T3                 = T_air
% 2) Router*Q - T3 + T2        = 0
% 3) R1*Q     - T2 + T1        = 0
% 4) Rliner*Q         - T1     = -T_LH2

Aeq = [ R2,      1,     0,    0;    % T_air = T3 + Q*R2
        Router, -1,     1,    0;    % T3 = T2 + Q*Router
        R1,      0,    -1,    1;    % T2 = T1 + Q*R1
        Rliner,  0,     0,   -1];   % T1 = T_LH2 + Q*Rliner

beq = [ T_air;
         0;
         0;
        -T_LH2 ];

x = Aeq \ beq;

Q  = x(1);
T3 = x(2);   % AIR-SIDE outer wall temp
T2 = x(3);   % Metal wall inner surface (outer tank inner wall)
T1 = x(4);   % Liner outer wall

Q_total = 1.13 * Q * 1.3;

% Equivalent h based on liner area
A_outer_liner = 1015456.009;    % mm^2
h_eq = Q_total / ((T2 - T1) * A_outer_liner);

% Results
fprintf('Q_total = %.6f W\n', Q_total);
fprintf('T3 = %.6f K (air-side outer wall)\n', T3);
fprintf('T2 = %.6f K (outer wall inner surface)\n', T2);
fprintf('T1 = %.6f K (liner outer surface)\n', T1);
fprintf('h_eq = %.6e W/mm^2K\n', h_eq);


%% calculate boil off

function compute_heat_flux(csvfile)

    % Load input geometry table
    T = readtable(csvfile);

    % Extract columns
    a_vals  = T.a_mm;
    b_vals  = T.b_mm;
    c_vals  = T.c_mm;
    n1_vals = T.n1;
    n2_vals = T.n2;

    % Constants (units mm, W/mm-K)
    t_ins   = 16;      
    t_outer = 2;
    k_liner = 0.0306;
    k_outer = 0.0306;
    k_ins   = 3.03e-08;
    conv_coeff = 10 / 1e6;      % W/mm^2-K
    T_air = 300;
    T_LH2 = 20;

    % Storage
    Q_totals = zeros(height(T),1);
    h_effs   = zeros(height(T),1);

    fprintf("Processing %d geometries...\n", height(T));

    for i = 1:height(T)

        a = a_vals(i);
        b = b_vals(i);
        c = c_vals(i);
        n1 = n1_vals(i);
        n2 = n2_vals(i);

        % Areas
        A_liner = superellipsoid_area(a, b, c, n1, n2);
        A_ins   = superellipsoid_area(a+t_ins, b+t_ins, c+t_ins, n1, n2);
        A_outer = superellipsoid_area(a+t_ins+t_outer, b+t_ins+t_outer, ...
                                      c+t_ins+t_outer, n1, n2);

        % Volumes for shape factor
        V_liner_inner = superellipsoid_volume(a, b, c, n1, n2);
        V_liner_outer = superellipsoid_volume(a+t_ins, b+t_ins, c+t_ins, n1, n2);
        V_outer_outer = superellipsoid_volume(a+t_ins+t_outer, b+t_ins+t_outer, ...
                                              c+t_ins+t_outer, n1, n2);

        V_ins   = V_liner_outer - V_liner_inner;
        V_outer = V_outer_outer - V_liner_outer;

        % Shape factors
        S_liner = shape_factor(A_liner, V_liner_outer - V_liner_inner, t_ins, a, b, c);
        S_ins   = shape_factor(A_ins,   V_ins,  t_ins,   a, b, c);
        S_outer = shape_factor(A_outer, V_outer, t_outer, a, b, c);

        % Compute heat flux & effective h
        [Q_total, T3, T2, T1, h_eq] = equivalent_heat_coeff( ...
            T_air, T_LH2, ...
            A_outer, ...
            conv_coeff, ...
            S_ins, k_ins, ...
            S_liner, k_liner, ...
            S_outer, k_outer, ...
            A_ins );

        Q_totals(i) = Q_total;
        h_effs(i)   = h_eq;

        fprintf("%4d / %4d: Q_total = %.6f W, h_eq = %.6e W/mm^2-K\n", ...
                i, height(T), Q_total, h_eq);
    end

    T.Q_total_W      = Q_totals;
    T.h_equivalent   = h_effs;

    writetable(T, "heat_flux_results.csv");
    fprintf("Saved output → heat_flux_results.csv\n");

end

function A = superellipsoid_area(a, b, c, n1, n2, N)
    if nargin < 6, N = 5250; end     % use smaller grid for speed

    phi   = linspace(-pi/2, pi/2, N);
    theta = linspace(-pi,   pi,   N);
    [PHI, THETA] = meshgrid(phi, theta);

    spow = @(x,p) sign(x) .* abs(x).^p;

    X = a * spow(cos(PHI), n1) .* spow(cos(THETA), n2);
    Y = b * spow(cos(PHI), n1) .* spow(sin(THETA), n2);
    Z = c * spow(sin(PHI), n1);

    [dX_dphi,   dX_dtheta]   = gradient(X, phi,   theta);
    [dY_dphi,   dY_dtheta]   = gradient(Y, phi,   theta);
    [dZ_dphi,   dZ_dtheta]   = gradient(Z, phi,   theta);

    r_phi   = cat(3, dX_dphi,   dY_dphi,   dZ_dphi);
    r_theta = cat(3, dX_dtheta, dY_dtheta, dZ_dtheta);

    cross_prod = cross(r_phi, r_theta, 3);
    dA = sqrt(sum(cross_prod.^2, 3));

    dphi   = phi(2)   - phi(1);
    dtheta = theta(2) - theta(1);

    A = sum(dA(:)) * dphi * dtheta;
end

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

function S = shape_factor(A, V_inner, thick, a, b, c)

    S_inf = 3.51 * sqrt(A);
    S_0   = A / thick;
    ls    = 2 * max([a, b, c]);

    expr = 1.26 - (2 - sqrt(A)/ls) / (9 * sqrt(1 - 4.79*(V_inner^(2/3)) / A));
    n = max(expr, 1.0);

    S = (S_0^n + S_inf^n)^(1/n);
end

function [Q_total, T3, T2, T1, h_eq] = equivalent_heat_coeff( ...
    T_air, T_LH2, ...
    A_outer, hc_mm, ...
    S_ins, k_ins, ...
    S_liner, k_liner, ...
    S_outer, k_outer, ...
    A_outer_liner)

    R2      = 1 / (A_outer * hc_mm);     % convection
    R1      = 1 / (k_ins   * S_ins);     % insulation
    Rliner  = 1 / (k_liner * S_liner);
    Router  = 1 / (k_outer * S_outer);

    Aeq = [
        R2,      1,     0,    0;
        Router, -1,     1,    0;
        R1,      0,    -1,    1;
        Rliner,  0,     0,   -1
    ];

    beq = [T_air; 0; 0; -T_LH2];

    x = Aeq \ beq;           % Solve linear system
    Q = x(1); T3 = x(2); T2 = x(3); T1 = x(4);

    Q_total = 1.13 * Q * 1.3;    % LH2 correction
    h_eq = Q_total / ((T_air - T_LH2) * A_outer_liner);
end

compute_heat_flux("superellipsoid_parametric_162.csv");

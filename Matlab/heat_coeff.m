% Known temperatures
T_air = 300;     % K
T_LH2 = 20;      % K

%% OUTER CONVECTION
A  = 1176179.901067;       % mm^2 (outer surface)
hc = 10;             % W/m^2K
hc_mm = hc/1e6;      % W/mm^2-K
R2 = 1 / (A * hc_mm);

%% INSULATION
k_eff_mW = 0.0303;        % mW/mmK
k_eff = k_eff_mW * 1e-6;  % W/mmK
S_ins = 67003.022;  % mm
G1 = k_eff * S_ins;
R1 = 1 / G1;

%% LINER
k_al = 0.0306;        % W/mmK
S_liner = 502540.155;  % mm
G_liner = k_al * S_liner;
Rliner = 1 / G_liner;

%% OUTER WALL METAL RESISTANCE
k_outer = 0.0306;       % W/mmK
S_outer = 582749.652;       % mm
Gouter  = k_outer * S_outer;
Router  = 1 / Gouter;

%% Unknowns vector: x = [Q; T3; T2; T1]
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

% Known temps
T_air = 293;    % K
T_LH2 = 20;     % K

% Outer film
A  = 809702.77;           % mm^2
hc = 10;                  % W/m^2K
hc_mm = hc/1e6;           % W/mm^2K
R2 = 1/(A * hc_mm);       % K/W

% Insulation
k_eff_mW = 0.3639;        % mW/mK
k_eff = k_eff_mW * 1e-6;  % W/mmK
S_ins = 147369.65;        % mm
G1 = k_eff * S_ins;       % W/K
R1 = 1 / G1;              % K/W

% Liner
k_al = 0.0306;            % W/mmK
S_liner = 138749.34;      % mm
G_liner = k_al * S_liner; % W/K
Rliner = 1 / G_liner;     % K/W

% Unknowns vector: x = [Q; T2; T1]
Aeq = [ R2,     1,    0;    % R2*Q + T2           = T_air
        R1,    -1,    1;    % R1*Q - T2 + T1      = 0
        Rliner, 0,   -1 ];  % Rliner*Q - T1       = -T_LH2
beq = [ T_air;  0;   -T_LH2 ];

x = Aeq \ beq;
Q  = x(1);
T2 = x(2);
T1 = x(3);

Q_LH2 = 1.13 * Q * 1.3;
T2_LH2 = T_air - Q_LH2*R2;
T1p_LH2 = 20 + Q_LH2*Rliner;

fprintf('Q   = %.6f W\n', Q_LH2);
fprintf('T2  = %.6f K (outer tank - inner wall)\n', T2_LH2);
fprintf('T1 = %.6f K (liner outer surface)\n', T1p_LH2);

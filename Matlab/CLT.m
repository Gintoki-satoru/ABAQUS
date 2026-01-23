angles_deg = [0,45,-45,0, 0, -45, 45, 0];
Nply = numel(angles_deg);
R = 241.09;     % Radius [mm]
t_tot     = 0.16*Nply;     % total thickness [mm]
R_mid = 241.09 + t_tot/2;
Pi    = 1;      
Po    = 0.0;

% Ply thicknesses [mm]
t_ply = t_tot / Nply;
tply  = t_ply * ones(1,Nply);

% UD lamina properties (plane stress) in MPa
E1  = 171000;    % [MPa]
E2  = 8530;     % [MPa]
G12 = 5630;      % [MPa]
nu12 = 0.287;
% (nu21 derived)
nu21 = nu12*E2/E1;

% ----------------------- MEMBRANE LOADS -----------------------
p = (Pi - Po);                    % [MPa]
% sigma_phi   = p*R_mid/(2*t_tot);  % [MPa]
% sigma_theta = sigma_phi;          % [MPa]
% tau_phitheta = 0.0;
% 
% % Membrane resultants [MPa*mm = N/mm]
% Nx  = sigma_phi   * t_tot;
% Ny  = sigma_theta * t_tot;
% Nxy = tau_phitheta * t_tot;

% Nvec = [Nx; Ny; Nxy];
% 
% % No moments in pure membrane loading
% Mvec = [0; 0; 0];

% Resultants from spherical equilibrium (valid regardless of laminate)
Nphi   = p*R_mid/2;          % MPa*mm = N/mm
Ntheta = p*R_mid/2;          % MPa*mm
Nxy    = 0;

Nvec = [Nphi; Ntheta; Nxy];
Mvec = [0; 0; 0];

% ----------------------- BUILD ABD -----------------------
Q = Qbar_UD_plane_stress(E1,E2,G12,nu12,nu21);

z = zeros(1,Nply+1);
z(1) = -t_tot/2;
for k=1:Nply
    z(k+1) = z(k) + tply(k);
end

A = zeros(3); B = zeros(3); D = zeros(3);
Qbar_all = zeros(3,3,Nply);

for k=1:Nply
    th = deg2rad(angles_deg(k));
    Qbar = Qbar_from_Q(Q, th);
    Qbar_all(:,:,k) = Qbar;
    zk  = z(k);
    zk1 = z(k+1);
    A = A + Qbar*(zk1 - zk);
    B = B + 0.5*Qbar*(zk1^2 - zk^2);
    D = D + (1/3)*Qbar*(zk1^3 - zk^3);
end

ABD = [A B; B D];

% ----------------------- SOLVE FOR MID-PLANE STRAINS/CURVATURES -----------------------
% [N; M] = [A B; B D]*[eps0; kappa]
NM = [Nvec; Mvec];
x = ABD \ NM;
eps0   = x(1:3);   % mid-plane strains [ - ]
kappa  = x(4:6);   % curvatures [1/mm]

% ----------------------- PLY STRESSES THROUGH THICKNESS -----------------------
% Evaluate at ply midsurfaces
z_mid = 0.5*(z(1:end-1) + z(2:end));
sigma12 = zeros(3,Nply);   % [sigma1; sigma2; tau12] in lamina coords
sigmaxy = zeros(3,Nply);   % [sigmax; sigmay; tauxy] in global x-y

for k=1:Nply
    th = deg2rad(angles_deg(k));
    Qbar = Qbar_all(:,:,k);
    % Global strain at this z
    eps_xy = eps0 + z_mid(k)*kappa;           % [ex; ey; gxy]
    % Global stress (x-y)
    sig_xy = Qbar * eps_xy;                   % [sx; sy; txy]
    sigmaxy(:,k) = sig_xy;

    % Convert global stress to lamina (1-2) stresses
    sigma12(:,k) = stress_xy_to_12(sig_xy, th);
end

% ----------------------- REPORT -----------------------
fprintf('--- CLT results (sphere membrane loading) ---\n');
fprintf('sigma_phi = sigma_theta = %.3f MPa (thin-sphere membrane)\n', sigma_phi);
fprintf('Nx = %.3f MPa*mm, Ny = %.3f MPa*mm\n', Nx, Ny);
fprintf('Mid-plane strains eps0 = [%.6e  %.6e  %.6e]\n', eps0(1), eps0(2), eps0(3));
fprintf('Curvatures kappa = [%.6e  %.6e  %.6e] 1/mm\n', kappa(1), kappa(2), kappa(3));

% Print ply lamina stresses
fprintf('\nPly mid-surface stresses in lamina coords (1-2):\n');
fprintf('Ply  Angle   sigma1(MPa)   sigma2(MPa)   tau12(MPa)\n');
for k=1:Nply
    fprintf('%3d  %6.1f  %12.3f  %12.3f  %12.3f\n', k, angles_deg(k), ...
        sigma12(1,k), sigma12(2,k), sigma12(3,k));
end

% ----------------------- PLOTS -----------------------
figure; hold on; grid on;
plot(z_mid, sigmaxy(1,:), '-o', 'LineWidth', 1.5);
plot(z_mid, sigmaxy(2,:), '-o', 'LineWidth', 1.5);
plot(z_mid, sigmaxy(3,:), '-o', 'LineWidth', 1.5);
xlabel('z (mm)'); ylabel('Stress (MPa)');
title('Global (x-y) ply mid-surface stresses from CLT');
legend('\sigma_x (meridional)','\sigma_y (hoop)','\tau_{xy}','Location','best');

figure; hold on; grid on;
plot(z_mid, sigma12(1,:), '-o', 'LineWidth', 1.5);
plot(z_mid, sigma12(2,:), '-o', 'LineWidth', 1.5);
plot(z_mid, sigma12(3,:), '-o', 'LineWidth', 1.5);
xlabel('z (mm)'); ylabel('Stress (MPa)');
title('Lamina (1-2) ply mid-surface stresses from CLT');
legend('\sigma_1','\sigma_2','\tau_{12}','Location','best');

%% ======================= FUNCTIONS =======================
function Q = Qbar_UD_plane_stress(E1,E2,G12,nu12,nu21)
% Plane-stress reduced stiffness Q in lamina coordinates (1-2)
den = 1 - nu12*nu21;
Q11 = E1/den;
Q22 = E2/den;
Q12 = nu12*E2/den;
Q66 = G12;
Q = [Q11 Q12 0;
     Q12 Q22 0;
     0   0   Q66];
end

function Qbar = Qbar_from_Q(Q, th)
% Transformed reduced stiffness Qbar for angle th (radians)
m = cos(th); n = sin(th);
Q11=Q(1,1); Q22=Q(2,2); Q12=Q(1,2); Q66=Q(3,3);

m2=m*m; n2=n*n; m3=m2*m; n3=n2*n; m4=m2*m2; n4=n2*n2;

Qbar11 = Q11*m4 + 2*(Q12+2*Q66)*m2*n2 + Q22*n4;
Qbar22 = Q11*n4 + 2*(Q12+2*Q66)*m2*n2 + Q22*m4;
Qbar12 = (Q11+Q22-4*Q66)*m2*n2 + Q12*(m4+n4);
Qbar16 = (Q11 - Q12 - 2*Q66)*m3*n - (Q22 - Q12 - 2*Q66)*m*n3;
Qbar26 = (Q11 - Q12 - 2*Q66)*m*n3 - (Q22 - Q12 - 2*Q66)*m3*n;
Qbar66 = (Q11 + Q22 - 2*Q12 - 2*Q66)*m2*n2 + Q66*(m4+n4);

Qbar = [Qbar11 Qbar12 Qbar16;
        Qbar12 Qbar22 Qbar26;
        Qbar16 Qbar26 Qbar66];
end

function sig12 = stress_xy_to_12(sigxy, th)
% Transform global stresses [sx; sy; txy] to lamina stresses [s1; s2; t12]
m = cos(th); n = sin(th);
sx = sigxy(1); sy = sigxy(2); txy = sigxy(3);

s1  = m^2*sx + n^2*sy + 2*m*n*txy;
s2  = n^2*sx + m^2*sy - 2*m*n*txy;
t12 = -m*n*sx + m*n*sy + (m^2 - n^2)*txy;

sig12 = [s1; s2; t12];
end

%%
% sigmaxy : 3 x Nplies  (global stresses at ply mids: [sx; sy; txy]) in MPa
% z_mid   : 1 x Nplies  (ply mid-thickness locations, mm)
% theta   : 1 x Nplies  (ply angles in degrees, material 1-axis w.r.t global x)
%
Xt = 3179.2;
Xc = -1705.3;
Yt = 55.701;
Yc = -367.44;
S  = 199.11;
theta = angles_deg;
Nplies = size(sigmaxy,2);

% --- Tsai-Wu coefficients (plane stress 1-2) ---
H1  = (1/Xt) - (1/Xc);
H2  = (1/Yt) - (1/Yc);
H11 = 1/(Xt*Xc);
H22 = 1/(Yt*Yc);
H66 = 1/(S*S);
H12 = -0.5 * sqrt(H11*H22);

sig12   = zeros(3,Nplies);   % [s1; s2; t12]
TW_index = zeros(1,Nplies);

for k = 1:Nplies
    sx  = sigmaxy(1,k);
    sy  = sigmaxy(2,k);
    txy = sigmaxy(3,k);

    m = cosd(theta(k));
    n = sind(theta(k));

    % --- Stress transform: global (x,y) -> local (1,2) ---
    s1  = m^2*sx + n^2*sy + 2*m*n*txy;
    s2  = n^2*sx + m^2*sy - 2*m*n*txy;
    t12 = -m*n*sx + m*n*sy + (m^2 - n^2)*txy;

    sig12(:,k) = [s1; s2; t12];

    % --- Tsai-Wu failure index ---
    TW_index(k) = H1*s1 + H2*s2 + H11*s1^2 + H22*s2^2 + H66*t12^2 + 2*H12*s1*s2;
end

% --- Plot Tsai-Wu index ---
figure; hold on; grid on;
plot(z_mid, TW_index, '-o', 'LineWidth', 1.5);
yline(1.0, '--', 'LineWidth', 1.5);
xlabel('z (mm)'); ylabel('Tsai-Wu index');
title('Tsai-Wu failure index per ply (mid-surface)');
legend('Tsai-Wu', 'Failure = 1', 'Location', 'best');

% Optional: highlight failed plies
failed = find(TW_index >= 1.0);
if ~isempty(failed)
    fprintf('Failed plies (TW>=1): %s\n', mat2str(failed));
    fprintf('Max Tsai-Wu = %.3f at ply %d (z_mid = %.3f mm)\n', max(TW_index), failed(TW_index(failed)==max(TW_index)), z_mid(TW_index==max(TW_index)));
else
    fprintf('No ply failed (TW<1). Max Tsai-Wu = %.3f\n', max(TW_index));
end

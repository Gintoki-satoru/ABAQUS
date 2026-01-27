tply   = 0.16;
angles = [90,30,-30,90];
Nply = numel(angles_deg);
R_in = 241.09;     % Radius [mm]
t_tot     = tply*Nply;     % total thickness [mm]
R = R_in + t_tot/2;
p_i = 1;
K   = 1;     % head flattening factor used in Eq.(67)

% UD lamina properties in MPa - 20K car epx
E1  = 173930;    % [MPa]
E2  = 15690;     % [MPa]
G12 = 11244;      % [MPa]
nu12 = 0.433;
nu21 = nu12*E2/E1;

% E1  = 161000;    % [MPa]
% E2  = 11380;     % [MPa]
% G12 = 5200;      % [MPa]
% nu12 = 0.32;

% Axial coordinate to check results results
xmax = 240;       % [mm]
Nx   = 241;
x    = linspace(0, xmax, Nx).';   % x=0 is the junction/edge where Q0x,M0x applied

[A, D, z_int] = clt_AD(E1,E2,G12,nu12, angles, tply);

A11=A(1,1); A12=A(1,2); A22=A(2,2);
D11=D(1,1); D12=D(1,2); D22=D(2,2);

% (Eq.30)
beta_c = ((A22 - (A12^2)/A11) / (4*R^2*D11))^(1/4);

% (Eq.34)
w_p = (p_i*R^2) / (A11*A22 - A12^2) * (A12^2 - A11);

% -------------------- DISCONTINUITY LOADS Q0x, M0x ----------------------
% (Eq.67)
Q0x = -(p_i*K^2/8) * ( (4*R^2*A11*D11)/(A11*A22 - A12^2) )^(1/4);
M0x = 0.0;

% -------------------- CLOSED-FORM FIELDS (Eqs.46-50) --------------------
bcx = beta_c*x;
expTerm = exp(-bcx);

% Eq.(46)
w = w_p + expTerm./(2*beta_c^3*D11) .* ...
    ( beta_c*M0x*(sin(bcx) - cos(bcx)) - Q0x*cos(bcx) );

% (Eq.22)
Nx_res = p_i*R/2;

% Eq.(48)
Ntheta = p_i*R + (A22 - A12^2/A11) .* expTerm./(2*beta_c^3*R*D11) .* ...
    ( beta_c*M0x*(cos(bcx) - sin(bcx)) + Q0x*cos(bcx) );

% Eq.(49)
Mx = expTerm./beta_c .* ( beta_c*M0x*(cos(bcx)+sin(bcx)) + Q0x*sin(bcx) );
Mtheta = (D12/D11) .* Mx;

% Eq.(50)
Qx = expTerm .* ( Q0x*(cos(bcx) - sin(bcx)) - 2*beta_c*M0x*sin(bcx) );

% -------------------- CONVERT (N,M) -> ply stresses (CLT) ---------------
Nxy = zeros(size(x));
Mxy = zeros(size(x));

eps0 = (A \ [Nx_res*ones(size(x)).'; Ntheta.'; Nxy.']).';
kap  = (D \ [Mx.'; Mtheta.'; Mxy.']).';

nply = numel(angles);
z_mid = 0.5*(z_int(1:end-1) + z_int(2:end));

sigx_ply   = zeros(Nx, nply);
sigt_ply   = zeros(Nx, nply);
taux_ply   = zeros(Nx, nply);

for k = 1:nply
    th = angles(k);
    Qbar = Qbar_plane_stress(E1,E2,G12,nu12, th);

    z = z_mid(k);
    eps_z = eps0 + z*kap;      

    sig = (Qbar * eps_z.').';
    sigx_ply(:,k) = sig(:,1);
    sigt_ply(:,k) = sig(:,2);
    taux_ply(:,k) = sig(:,3);
end

% -------------------- PLOTS ---------------------------------------------
figure; plot(x, Ntheta, 'LineWidth', 1.5); grid on;
xlabel('x from junction [mm]'); ylabel('N_\theta [MPa·mm]');
title('Circumferential stress resultant N_\theta(x) (BTCS)');

figure; plot(x, Mx, 'LineWidth', 1.5); grid on;
xlabel('x from junction [mm]'); ylabel('M_x [MPa·mm^2]');
title('Axial bending moment M_x(x) (BTCS)');

figure; plot(x, Qx, 'LineWidth', 1.5); grid on;
xlabel('x from junction [mm]'); ylabel('Q_x [MPa·mm]');
title('Shear force Q_x(x) (BTCS)');

% Ply mid-surface hoop stress along x
figure; plot(x, sigt_ply(:,end), 'LineWidth', 1.5); grid on;
xlabel('x from junction [mm]'); ylabel('\sigma_\theta [MPa]');
title('Hoop stress at ply mid-surface (outermost ply)');

% Through-thickness ply mid-surface stresses at x=0
[~,ix0] = min(abs(x-0));
figure; plot(z_mid, sigt_ply(ix0,:), '-o', 'LineWidth', 1.5); grid on;
xlabel('z (thickness coord) [mm]'); ylabel('\sigma_\theta [MPa]');
title('Hoop stress through thickness at x=0 (ply mid-surfaces)');

%% ========================== FUNCTIONS ===================================
function [A, D, z_int] = clt_AD(E1,E2,G12,nu12, angles_deg, tply)
    % Plane-stress CLT A and D for a laminate
    n = numel(angles_deg);
    h = n*tply;
    z_int = linspace(-h/2, +h/2, n+1);  % interfaces

    A = zeros(3,3);
    D = zeros(3,3);

    for k = 1:n
        th = angles_deg(k);
        Qb = Qbar_plane_stress(E1,E2,G12,nu12, th);
        z1 = z_int(k);
        z2 = z_int(k+1);

        A = A + Qb*(z2 - z1);
        D = D + (1/3)*Qb*(z2^3 - z1^3);
    end
end

function Qb = Qbar_plane_stress(E1,E2,G12,nu12, theta_deg)
    % Orthotropic lamina plane-stress Qbar (global x-theta axes)
    nu21 = nu12*E2/E1;
    den  = 1 - nu12*nu21;

    Q11 = E1/den;
    Q22 = E2/den;
    Q12 = nu12*E2/den;
    Q66 = G12;

    m = cosd(theta_deg);
    n = sind(theta_deg);

    m2 = m*m; n2 = n*n;
    m3 = m2*m; n3 = n2*n;
    m4 = m2*m2; n4 = n2*n2;

    Qb11 = Q11*m4 + 2*(Q12+2*Q66)*m2*n2 + Q22*n4;
    Qb22 = Q11*n4 + 2*(Q12+2*Q66)*m2*n2 + Q22*m4;
    Qb12 = (Q11+Q22-4*Q66)*m2*n2 + Q12*(m4+n4);
    Qb16 = (Q11 - Q12 - 2*Q66)*m3*n - (Q22 - Q12 - 2*Q66)*m*n3;
    Qb26 = (Q11 - Q12 - 2*Q66)*m*n3 - (Q22 - Q12 - 2*Q66)*m3*n;
    Qb66 = (Q11+Q22-2*Q12-2*Q66)*m2*n2 + Q66*(m4+n4);

    Qb = [Qb11 Qb12 Qb16;
          Qb12 Qb22 Qb26;
          Qb16 Qb26 Qb66];
end

% Geometry
R = 180;     % Radius [mm]
t     = 28;     % total thickness [mm]
R_mid = 180 + t/2;

% Pressures
Pi = 70;         % internal pressure [MPa]
Po = 0.0;          % external pressure [MPa]

% Material properties - car_epx
% E1  = 171000;       
% E2  = 8530;
% E3  = 8530;
% G12 = 5630;
% G13 = 5630;
% G23 = 2640;
% nu12 = 0.287;
% nu13 = 0.287;
% nu23 = 0.369;

% Im7
% E1  = 161000;       
% E2  = 11380;
% E3  = 11380;
% G12 = 5200;
% G13 = 5200;
% G23 = 3900;
% nu12 = 0.32;
% nu13 = 0.32;
% nu23 = 0.45;

E1  = 10000;   % MPa
E2  = 147000;
E3  = 10000;
G12 = 7000;
G13 = 3700;
G23 = 7000;
nu21 = 0.27;
nu13 = 0.35;
nu23 = 0.27;
nu12 = nu21*E1/E2;

% Quasi-isotropic symmetric layup
angles_deg = [45 -45 45 -45 45 -45 45 -45 45 -45];
Nply = numel(angles_deg);
tply = t / Nply;

% Through-thickness evaluation resolution per ply
nPtsPerPly = 40;

% ---------------------- COMPUTE LOCAL STIFFNESS -------------------------
C_local = orthotropic_stiffness_3D(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23);

% Extract needed local terms
C11 = C_local(1,1);  C12 = C_local(1,2);  C13 = C_local(1,3);
C22 = C_local(2,2);  C23 = C_local(2,3);  C33 = C_local(3,3);

% ---------------------- COMPUTE K1, K2 -----------------
beta = -(C22 + 2*C23 + C33 - C12 - C13) / C11;
disc = 1 - 4*beta;
K1 = (-1 - sqrt(disc))/2;
K2 = (-1 + sqrt(disc))/2;

% ---------------------- SOLVE A1, A2 FROM RADIAL BCs --------------------
Hi = -t/2;   Ho = +t/2;     % thickness coordinate h
Ri = R_mid + Hi;            % inner radius
Ro = R_mid + Ho;            % outer radius

% For uniform expansion solution:
Q1 = (C12 + C13 + K1*C11);
Q2 = (C12 + C13 + K2*C11);

A = [ Q1*Ri^(K1-1),  Q2*Ri^(K2-1);
      Q1*Ro^(K1-1),  Q2*Ro^(K2-1) ];

b = [ -Pi; -Po ];
x = A \ b;
A1 = x(1); A2 = x(2);

% ---------------------- THROUGH-THICKNESS GRID BY PLY -------------------
h_all = [];
ply_id = [];
alpha_all = [];

for k = 1:Nply
    h_k0 = Hi + (k-1)*tply;
    h_k1 = Hi + (k)*tply;
    hk = linspace(h_k0, h_k1, nPtsPerPly).';
    h_all   = [h_all; hk];
    ply_id  = [ply_id; k*ones(size(hk))];
    alpha_all = [alpha_all; angles_deg(k)*ones(size(hk))];
end

h1_all = h_all + R_mid;

% ---------------------- COMPUTE STRESSES PLY-WISE -----------------------
sigma_phi      = zeros(size(h_all));
sigma_theta    = zeros(size(h_all));
sigma_h        = zeros(size(h_all));
tau_phitheta   = zeros(size(h_all));

for k = 1:Nply
    ak = deg2rad(angles_deg(k));
    Ck = rotate_C_about_1(C_local, ak);

    % transformed constants for that ply
    C11k = Ck(1,1); C12k = Ck(1,2); C13k = Ck(1,3);
    C22k = Ck(2,2); C23k = Ck(2,3); C33k = Ck(3,3);
    C16k = Ck(1,6); C26k = Ck(2,6); C36k = Ck(3,6);
    C14k = Ck(1,4);
    C24k = Ck(2,4);
    C34k = Ck(3,4);

    idx = (ply_id == k);
    h1  = h1_all(idx);

    p1 = h1.^(K1-1);
    p2 = h1.^(K2-1);

    sigma_h(idx)     = (C12k + C13k + K1*C11k)*A1.*p1 + (C12k + C13k + K2*C11k)*A2.*p2; % radial -> Abaqus S11
    sigma_phi(idx)   = (C22k + C23k + K1*C12k)*A1.*p1 + (C22k + C23k + K2*C12k)*A2.*p2; % meridional -> Abaqus S22
    sigma_theta(idx) = (C23k + C33k + K1*C13k)*A1.*p1 + (C23k + C33k + K2*C13k)*A2.*p2; % hoop -> Abaqus S33
    tau_phitheta(idx) = (K1*C14k + C24k + C34k)*A1.*p1 + (K2*C14k + C24k + C34k)*A2.*p2;
end

% ---------------------- PLOT RESULTS ------------------------------------
r_all = R_mid + h_all;

% Meridional stress
figure;
plot(h_all, sigma_phi, 'LineWidth', 1.5); grid on;
xlabel('Thickness coordinate h [mm]');
ylabel('\sigma_\phi [MPa]');
title('Meridional stress through thickness');
xline(-t/2,'--'); xline(t/2,'--');

% Hoop stress
figure;
plot(h_all, sigma_theta, 'LineWidth', 1.5); grid on;
xlabel('Thickness coordinate h [mm]');
ylabel('\sigma_\theta [MPa]');
title('Hoop stress through thickness');
xline(-t/2,'--'); xline(t/2,'--');

% Radial stress
figure;
plot(h_all, sigma_h, 'LineWidth', 1.5); grid on;
xlabel('Thickness coordinate h [mm]');
ylabel('\sigma_h [MPa]');
title('Radial stress through thickness');
xline(-t/2,'--'); xline(t/2,'--');

% In-plane shear stress
figure;
plot(h_all, tau_phitheta, 'LineWidth', 1.5); grid on;
xlabel('Thickness coordinate h [mm]');
ylabel('\tau_{\phi\theta} [MPa]');
title('In-plane shear stress through thickness');
xline(-t/2,'--'); xline(t/2,'--');

%% ---------------------- OPTIONAL: EXPORT TABLE --------------------------
T = table(r_all, h_all, ply_id, alpha_all, sigma_phi, sigma_theta, sigma_h, tau_phitheta);
disp(T(1:10,:));

%% ====================== FUNCTIONS =======================================

function C = orthotropic_stiffness_3D(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)
    % Builds 6x6 stiffness matrix in local material coordinates (1,2,3)
    % using compliance inversion.
    nu21 = nu12 * E2/E1;
    nu31 = nu13 * E3/E1;
    nu32 = nu23 * E3/E2;

    S = zeros(6,6);
    S(1,1) = 1/E1;    S(2,2) = 1/E2;    S(3,3) = 1/E3;
    S(1,2) = -nu21/E2; S(2,1) = -nu12/E1;
    S(1,3) = -nu31/E3; S(3,1) = -nu13/E1;
    S(2,3) = -nu32/E3; S(3,2) = -nu23/E2;

    S(4,4) = 1/G23;
    S(5,5) = 1/G13;
    S(6,6) = 1/G12;

    C = inv(S);
end

function Ck = rotate_C_about_1(C, alpha)
    % Rotate stiffness about axis-1 by angle alpha (engineering shear convention)
    m = cos(alpha); n = sin(alpha);

    % This rotates the (2-3) plane (meridional-hoop) about axis 1 (radial).
    T1 = [ 1, 0, 0, 0, 0, 0;
           0, m^2, n^2,  2*m*n, 0, 0;
           0, n^2, m^2, -2*m*n, 0, 0;
           0,-m*n, m*n,  (m^2-n^2),0,0;
           0, 0, 0, 0,  m, -n;
           0, 0, 0, 0,  n,  m ];

    T2 = [ 1, 0, 0, 0, 0, 0;
           0, m^2, n^2,  m*n, 0, 0;
           0, n^2, m^2, -m*n, 0, 0;
           0,-2*m*n,2*m*n,(m^2-n^2),0,0;
           0, 0, 0, 0,  m, -n;
           0, 0, 0, 0,  n,  m ];

    % Ck = T1(-a)*C*T2(a)
    nneg = -n; mneg = m;
    T1m = [ 1, 0, 0, 0, 0, 0;
            0, mneg^2, nneg^2,  2*mneg*nneg, 0, 0;
            0, nneg^2, mneg^2, -2*mneg*nneg, 0, 0;
            0,-mneg*nneg, mneg*nneg,(mneg^2-nneg^2),0,0;
            0, 0, 0, 0,  mneg, -nneg;
            0, 0, 0, 0,  nneg,  mneg ];

    Ck = T1m * C * T2;
end

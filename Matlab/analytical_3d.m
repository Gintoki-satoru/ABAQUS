data = readmatrix('3D.xlsx');
X  = data(:,1);  Y  = data(:,2);  Z  = data(:,3);
S11 = data(:,4); S22 = data(:,5); S33 = data(:,6);
S12 = data(:,7); S13 = data(:,8); S23 = data(:,9);

N = length(X);

% Radial coordinate
r = sqrt(X.^2 + Y.^2);

% ---- Sort along meridian (pick one: by Z or by r). For pole->equator, Z is good.
[Zs, idx] = sort(Z, 'descend');   % pole (high Z) to equator (lower Z)
X = X(idx); Y = Y(idx); Z = Z(idx);
S11=S11(idx); S22=S22(idx); S33=S33(idx);
S12=S12(idx); S13=S13(idx); S23=S23(idx);
r = r(idx);

% ---- Meridional tangent in r-z plane (use gradient for same-length output)
dr = gradient(r);
dz = gradient(Z);

tm_rz = [dr, dz];
tm_norm = sqrt(tm_rz(:,1).^2 + tm_rz(:,2).^2);
tm_rz = tm_rz ./ tm_norm;

% ---- Convert (r,z) tangent to 3D
eps_r = 1e-12;
r_safe = max(r, eps_r);

tm3D = [ tm_rz(:,1).*X./r_safe, ...
         tm_rz(:,1).*Y./r_safe, ...
         tm_rz(:,2) ];

tm3D = tm3D ./ vecnorm(tm3D,2,2);

% ---- Hoop direction (circumferential)
etheta = [ -Y./r_safe, X./r_safe, zeros(N,1) ];
etheta = etheta ./ vecnorm(etheta,2,2);

% ---- Compute stresses
sigma_m     = zeros(N,1);
sigma_hoop  = zeros(N,1);

for i = 1:N
    S = [ S11(i)  S12(i)  S13(i);
          S12(i)  S22(i)  S23(i);
          S13(i)  S23(i)  S33(i) ];

    m = tm3D(i,:)';
    h = etheta(i,:)';

    sigma_m(i)    = m' * S * m;
    sigma_hoop(i) = h' * S * h;
end

figure;
plot(r, sigma_m, 'LineWidth', 1.8);
grid on;
xlabel('\rho (mm)');
ylabel('\sigma_{meridional}');
title('Meridional Stress Distribution');

figure;
plot(r, sigma_hoop, 'LineWidth', 1.8);
grid on;
xlabel('\rho (mm)');
ylabel('\sigma_{hoop}');
title('Hoop Stress Distribution');
%% Sphere

% Inputs from your CSV/data matrix
X  = data(:,1);
Y  = data(:,2);
Z  = data(:,3);

S11 = data(:,4);  % sigma_xx
S22 = data(:,5);  % sigma_yy
S33 = data(:,6);  % sigma_zz
S12 = data(:,7);  % sigma_xy
S13 = data(:,8);  % sigma_xz
S23 = data(:,9);  % sigma_yz

% ---- If sphere center is not origin, uncomment and set center:
% Cx = 0; Cy = 0; Cz = 0;
% X = X - Cx; Y = Y - Cy; Z = Z - Cz;

% Spherical angles (physics convention)
r = sqrt(X.^2 + Y.^2 + Z.^2);
theta = acos(Z ./ r);         % polar angle [0..pi] from +Z
phi   = atan2(Y, X);          % azimuth [-pi..pi] from +X toward +Y

% Preallocate
N = numel(X);
S_rr  = zeros(N,1);
S_tt  = zeros(N,1);  % theta-theta (meridional)
S_pp  = zeros(N,1);  % phi-phi (hoop)
S_rt  = zeros(N,1);
S_rp  = zeros(N,1);
S_tp  = zeros(N,1);

for i = 1:N
    ct = cos(theta(i));  st = sin(theta(i));
    cp = cos(phi(i));    sp = sin(phi(i));

    % Unit vectors in Cartesian components
    e_r     = [ st*cp;  st*sp;  ct ];
    e_theta = [ ct*cp;  ct*sp; -st ];   % meridional direction
    e_phi   = [ -sp;     cp;    0  ];   % hoop direction

    % Global stress tensor at node i
    S = [ S11(i), S12(i), S13(i);
          S12(i), S22(i), S23(i);
          S13(i), S23(i), S33(i) ];

    % Transform: sigma_ab = e_a^T * S * e_b
    S_rr(i) = e_r.'     * S * e_r;
    S_tt(i) = e_theta.' * S * e_theta;  % meridional
    S_pp(i) = e_phi.'   * S * e_phi;    % hoop

    S_rt(i) = e_r.'     * S * e_theta;
    S_rp(i) = e_r.'     * S * e_phi;
    S_tp(i) = e_theta.' * S * e_phi;
end

% Outputs you likely want:
S_meridional = S_tt;   % sigma_theta-theta
S_hoop       = S_pp;   % sigma_phi-phi

r = sqrt(X.^2 + Y.^2 + Z.^2);

R_mid = mean(r);   % or use known value, e.g., 200

% Sort by theta
[theta_sorted, idx] = sort(theta);

% True arc length along meridian
s = R_mid * cos(theta_sorted);   % arc length

% Plot meridional stress
figure;
plot(s, S_meridional(idx), 'LineWidth', 1.5);
grid on;
xlabel('Arc length s (mm)');
ylabel('\sigma_{\theta\theta}');
title('Meridional stress vs true surface distance');

% Plot hoop stress
figure;
plot(s, S_hoop(idx), 'LineWidth', 1.5);
grid on;
xlabel('Arc length s (mm)');
ylabel('\sigma_{\phi\phi}');
title('Hoop stress vs true surface distance');
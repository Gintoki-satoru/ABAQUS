% Load data
data = readmatrix('3D.xlsx');

X  = data(:,1);
Y  = data(:,2);
Z  = data(:,3);

S11 = data(:,4);
S22 = data(:,5);
S33 = data(:,6);
S12 = data(:,7);
S13 = data(:,8);
S23 = data(:,9);

N = length(X);

% Radial coordinate
r = sqrt(X.^2 + Y.^2);

% --- Meridional direction (same as before) ---
dr = diff(r);
dz = diff(Z);

tm = [dr, dz];
tm = [tm; tm(end,:)];          % match length

tm_norm = sqrt(tm(:,1).^2 + tm(:,2).^2);
tm = tm ./ tm_norm;

% Convert (r,z) → (x,y,z)
tm3D = [ (tm(:,1).*X./r), ...
         (tm(:,1).*Y./r), ...
          tm(:,2) ];

for i = 1:N
    tm3D(i,:) = tm3D(i,:) / norm(tm3D(i,:));
end

% --- Circumferential (hoop) direction ---
etheta = [ -Y./r, X./r, zeros(N,1) ];

for i = 1:N
    etheta(i,:) = etheta(i,:) / norm(etheta(i,:));
end

% Build stress tensor + compute stresses
sigma_m     = zeros(N,1);
sigma_theta = zeros(N,1);

for i = 1:N
    
    % Stress tensor
    S = [ S11(i)  S12(i)  S13(i);
          S12(i)  S22(i)  S23(i);
          S13(i)  S23(i)  S33(i) ];
    
    tm_vec = tm3D(i,:)';
    th_vec = etheta(i,:)';
    
    % Meridional stress
    sigma_m(i) = tm_vec' * S * tm_vec;
    
    % Circumferential (hoop) stress
    sigma_theta(i) = th_vec' * S * th_vec;
end

% Arc length for plotting
ds = sqrt(dr.^2 + dz.^2);
s = [0; cumsum(ds)];

% --- Plot both ---
figure; hold on; grid on; box on;
plot(X, sigma_m, 'LineWidth', 2);
plot(X, sigma_theta, 'LineWidth', 2);

xlabel('Arc length s (mm)');
ylabel('Stress (MPa)');
title('Meridional and Circumferential Stress Along Ellipsoid Surface');
legend('\sigma_m (meridional)', '\sigma_\theta (circumferential)');

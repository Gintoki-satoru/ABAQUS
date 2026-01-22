% --- Abaqus extracted data: [h(mm)  S_rr  S_meridional  S_hoop  tau_merid_hoop] ---
abaqus_data = [
    0.000000   -0.0999638   17.2204    1.52995   -0.000157163
    0.000000   -0.1000500   17.2187    1.52979    2.63E-05
    0.159995   -0.0831221   17.1997    1.53472   -0.000158541
    0.159995   -0.0832509   17.1998   13.2846     6.78383
    0.159995   -0.0831082    5.45091  13.2839     6.78318
    0.159995   -0.0832021    5.45125  13.2669     6.77086
    0.319967   -0.0664452    5.44891  13.2694     6.77245
    0.319967   -0.0664298    5.44978  13.2658     6.77001
    0.319967   -0.0664409    5.44869  13.2686    -6.77189
    0.319967   -0.0664609    5.44964  13.2509    -6.75874
    0.479956   -0.0497537    5.44713  13.2503    -6.75839
    0.479956   -0.0497660    5.44693  13.2375    -6.74820
    0.479956   -0.0497727    5.44577  13.2377    -6.74833
    0.639940   -0.0331112    5.44591  13.2365    -6.74755
    0.639940   -0.0330839    5.44556  13.2367     6.74767
    0.639940   -0.0331241    5.44349  13.2219     6.73600
    0.639940   -0.0330803   5.44349   13.2219     6.73600
    0.799916   -0.0164875   5.44349   13.2219     6.73600
    0.799916   -0.0165444   17.0873   1.52995   -0.000214224
    0.799916   -0.0165850   17.0873   1.52995   -0.000214224
    0.799916   -0.0164516   17.0873   1.52995   -0.000214224
    0.959912    4.60E-05    17.0873   1.52995   -0.000214224
    0.959912    7.25E-05    17.0873   1.52995   -0.000214224
];

% --- Extract vectors ---
h_all_abaqus       = abaqus_data(:,1);
sigma_theta_abaq = abaqus_data(:,4);
sigma_merid_abaq = abaqus_data(:,3);
sigma_tau_abaq = abaqus_data(:,5);

mask = ~isnan(h_all_abaqus) & ~isnan(sigma_theta_abaq);
h_all_abaqus = h_all_abaqus(mask);
sigma_theta_abaq = sigma_theta_abaq(mask);
sigma_merid_abaq = sigma_merid_abaq(mask);
sigma_tau_abaq = sigma_tau_abaq(mask);

[h_all_abaqus, idx] = sort(h_all_abaqus);
sigma_theta_abaq  = sigma_theta_abaq(idx);
sigma_merid_abaq  = sigma_merid_abaq(idx);
sigma_tau_abaq  = sigma_tau_abaq(idx);

% Total thickness from Abaqus data
t = max(h_all_abaqus);

% Center Abaqus thickness coordinate about mid-plane
h_all_c = h_all_abaqus - t/2;

% hoop stress##########################################
figure; hold on; grid on;
% CLT hoop stress (ply mid-surfaces)
plot(z_mid, sigmaxy(2,:), '-o', ...
     'LineWidth', 1.5, ...
     'DisplayName', 'CLT (ply mid-surface)');

% Abaqus hoop stress (integration points)
plot(h_all_c, sigma_theta_abaq, '-','DisplayName', 'Abaqus');

% Marie sol
plot(h_all, sigma_theta,'--','DisplayName', 'Marie');

xlabel('Thickness coordinate z [mm]');
ylabel('\sigma_\theta [MPa]');
title('Hoop stress through thickness');
legend('Location','best');

% Merid stress##########################################
figure; hold on; grid on;
% CLT hoop stress (ply mid-surfaces)
plot(z_mid, sigmaxy(1,:), '-o', ...
     'LineWidth', 1.5, ...
     'DisplayName', 'CLT (ply mid-surface)');

% Abaqus hoop stress (integration points)
plot(h_all_c, sigma_merid_abaq, '-','DisplayName', 'Abaqus');

% Marie sol
plot(h_all, sigma_phi,'--','DisplayName', 'Marie');

xlabel('Thickness coordinate z [mm]');
ylabel('\sigma_\phi [MPa]');
title('Meridional stress through thickness');
legend('Location','best');

% tau12 stress##########################################
figure; hold on; grid on;
% CLT hoop stress (ply mid-surfaces)
plot(z_mid, sigmaxy(3,:), '-o', ...
     'LineWidth', 1.5, ...
     'DisplayName', 'CLT (ply mid-surface)');

% Abaqus hoop stress (integration points)
plot(h_all_c, sigma_tau_abaq, '-','DisplayName', 'Abaqus');

% Marie sol
plot(h_all, tau_phitheta,'--','DisplayName', 'Marie');

xlabel('Thickness coordinate z [mm]');
ylabel('\tau_{\phi\theta} [MPa]');
title('In-plane shear stress through thickness');
legend('Location','best');
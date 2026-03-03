% ---- Data from your table (3 angles) ----
phi = [15; 45; 75];

% Global (show one representative component ~188 MPa)
S33 = [188.2; 188.2949; 188.28];

% Local ply stresses
S11p = [218.08; 340.8925; 219.8881];
S33p = [156.6026;  35.6978; 156.8818];
S12p = [54.5863;   0.000229; 54.7575];

% Failure indices
TW   = [2.200091; 0.924933; 2.20219];
Puck_IFF = [2.826037; 0.640884; 2.831114];

% Allowables (MPa)
Xt = 3179.2;  Xc = -1705.3;
Yt = 55.701;  Yc = -367.44;
S  = 199.11;

% ---- Utilization ratios (dimensionless) ----
U1  = S11p / Xt;                % (all tensile here)
U2  = S33p / Yt;                % (all tensile here)
U12 = abs(S12p) / S;

% ---- Figure: 3 panels ----
figure;

% Panel A: global same
subplot(3,1,1);
plot(phi, S33, 'o-','LineWidth',1.5);
ylabel('Global S_{22} (MPa)');
grid on;
title('Same global stress, different local ply stress and failure index');

% Panel B: local components (or swap to U1/U2/U12 if you prefer)
subplot(3,1,2);
bar(phi, [S11p S33p S12p], 'grouped');
ylabel('Local ply stress (MPa)');
legend('S_{11}^p (fiber)','S_{22}^p (trans)','S_{12}^p (shear)','Location','best');
grid on;

% subplot(3,1,2);
% bar(phi, [U1 U2 U12], 'grouped'); hold on;
% yline(1.0,'--','LineWidth',1.2);
% ylabel('Utilization ratio');
% legend('\sigma_1/X_t','\sigma_2/Y_t','|\tau_{12}|/S','= 1','Location','best');
% grid on;

% Panel C: failure indices
subplot(3,1,3);
plot(phi, TW, 'o-','LineWidth',1.5); hold on;
plot(phi, Puck_IFF, 's-','LineWidth',1.5);
ylabel('Failure Index');
xlabel('\phi (deg)');
legend('Tsai-Wu','Puck IFF','FI = 1','Location','best');
grid on;

% ################################################################

% Utilization ratios
U1  = S11p ./ Xt;
U2  = S22p ./ Yt;
U12 = abs(S12p) ./ S;

% Data matrix
data = [S11p S22p S12p];

% -------- Plot --------
figure;
b = bar(phi, data, 'grouped');
hold on;
grid on;
ylim([-60 375])

xlabel('Fiber angle (deg)');
ylabel('Local ply stress (MPa)');

legend({'S_{11}^p (fiber)', ...
        'S_{22}^p (transverse)', ...
        '|S_{12}^p| (shear)'}, ...
        'Location','northoutside','Orientation','horizontal');

glob = 188;
yline(Yt,'--','Transverse limit - \frac{$S_{22}^p}{$Y_t$}$','Interpreter','latex','LineWidth',1.2, HandleVisibility='off');
yline(S ,'--','Shear limit - \frac{|$S_{12}^p|}{S}$','Interpreter','latex','LineWidth',1.2, HandleVisibility='off');
h = yline(glob, '.', 'Hoop stress = 188 MPa', ...
          'LineWidth',1.2, ...
          'HandleVisibility','off');

% Move text to right side
h.LabelHorizontalAlignment = 'right';
h.LabelVerticalAlignment   = 'bottom';
title('Local Ply Stresses and Corresponding Stress-to-Strength Ratios');

% -------- Add utilization values on top of bars --------
ngroups = size(data,1);
nbars   = size(data,2);

% Get bar positions
for i = 1:nbars
    x = b(i).XEndPoints;
    y = b(i).YEndPoints;
    
    if i == 1
        labels = U1;
    elseif i == 2
        labels = U2;
    else
        labels = U12;
    end
    
    for j = 1:ngroups
        text(x(j), y(j), sprintf('%.2f', labels(j)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',9);
    end
end
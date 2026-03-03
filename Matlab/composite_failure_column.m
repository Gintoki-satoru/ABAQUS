% ============================================================
% n = 0.65
% ============================================================
strain_limit = 0.005;

layups = {
    '(90/15/-15/90)'
    '(90/30/-30/90)'
    '(90/45/-45/90)'
    '(90/60/-60/90)'
    '(90/75/-75/90)'
    '(60/-30/-60/30)'
    '(15/-15/15/-15)'
    '(30/-30/30/-30)'
    '(45/-45/45/-45)'
    '(60/-60/60/-60)'
    '(75/-75/75/-75)'
};

TW = [0.6969; 0.8110; 0.9897; 1.1932; 1.3515; 0.7751; 2.0432; 1.5408; 0.8193; 1.1266; 1.3364];
PK = [0.4608; 0.5693; 0.7495; 1.0690; 1.3744; 0.5335; 2.5004; 1.4312; 0.5964; 0.8924; 1.2810];
eps = [0.00164;0.00200;0.00310;0.00480;0.00600;0.00190;0.01020;0.00540;0.00226;0.00340;0.00544];

eps_norm = eps / strain_limit;  % normalize only by strain limit

% ============================================================
% 1) FIRST PLOT: only the 5 angle-ply cases (90/±theta/90)
% ============================================================
idx1 = 1:6;

Y1 = [TW(idx1), PK(idx1), eps_norm(idx1)];
labels1 = layups(idx1);

figure('Color','w','Position',[200 120 1050 520]);
b1 = bar(Y1, 'grouped'); grid on; box on;
ylim([0, max([Y1(:); 1]) * 1.15]);

set(gca, 'XTick', 1:numel(idx1), 'XTickLabel', labels1, 'FontSize', 11);
xtickangle(20);

ylabel('Failure index / Normalized strain');
title(sprintf('Multidirectional Laminate for Super Ellipsoid with n_1 = 0.65'));

legend({'Tsai-Wu','Puck','Max Strain index'}, 'Location','northeast');

% failure/limit line at 1
hold on;
yline(1, '--', 'LineWidth', 1.2, HandleVisibility='off');
hold off;

% optional: value labels on bars
for k = 1:numel(b1)
    xt = b1(k).XEndPoints;
    yt = b1(k).YEndPoints;
    vals = b1(k).YData;
    text(xt, yt, compose('%.3f', vals), 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize', 9);
end

% ============================================================
% 2) SECOND PLOT: remaining 6 layups
% ============================================================
idx2 = 7:11;

Y2 = [TW(idx2), PK(idx2), eps_norm(idx2)];
labels2 = layups(idx2);

figure('Color','w','Position',[200 120 1050 520]);
b2 = bar(Y2, 'grouped'); grid on; box on;
ylim([0, max([Y2(:); 1]) * 1.15]);

set(gca, 'XTick', 1:numel(idx2), 'XTickLabel', labels2, 'FontSize', 11);
xtickangle(20);

ylabel('Failure index / Normalized strain');
title(sprintf('Angle ply Laminate for Super Ellipsoid with n_1 = 0.65'));

legend({'Tsai-Wu','Puck','Max Strain index'}, 'Location','northeast');

hold on;
yline(1, '--', 'LineWidth', 1.2, 'HandleVisibility','off');
hold off;

% optional: value labels on bars
for k = 1:numel(b2)
    xt = b2(k).XEndPoints;
    yt = b2(k).YEndPoints;
    vals = b2(k).YData;
    text(xt, yt, compose('%.3f', vals), 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize', 9);
end

%% ============================================================
% sphere
% ============================================================
strain_limit = 0.005;

layups = {
    '(90/15/-15/90)'
    '(90/30/-30/90)'
    '(90/45/-45/90)'
    '(90/60/-60/90)'
    '(90/75/-75/90)'
    '(60/-30/-60/30)'
    '(15/-15/15/-15)'
    '(30/-30/30/-30)'
    '(45/-45/45/-45)'
    '(60/-60/60/-60)'
    '(75/-75/75/-75)'
};

TW = [
    0.925     % (90/15/-15/90)
    1.14      % (90/30/-30/90)
    1.49      % (90/45/-45/90)
    1.9926    % (90/60/-60/90)
    2.2684    % (90/75/-75/90)
    1.0364    % (60/-30/-60/30)
    2.3146    % (15/-15/15/-15)
    1.9688    % (30/-30/30/-30)
    0.9302    % (45/-45/45/-45)
    1.6583    % (60/-60/60/-60)
    2.209     % (75/-75/75/-75)
];

PK = [
    0.6414    % (90/15/-15/90)
    0.953     % (90/30/-30/90)
    1.49      % (90/45/-45/90)
    2.2685    % (90/60/-60/90)
    3.1732    % (90/75/-75/90)
    0.7398    % (60/-30/-60/30)
    3.2044    % (15/-15/15/-15)
    1.9464    % (30/-30/30/-30)
    0.6511    % (45/-45/45/-45)
    1.5254    % (60/-60/60/-60)
    2.8689    % (75/-75/75/-75)
];

eps = [
    0.0021     % (90/15/-15/90)
    0.0030     % (90/30/-30/90)
    0.0060     % (90/45/-45/90)
    0.009571   % (90/60/-60/90)
    0.0130     % (90/75/-75/90)
    0.0024     % (60/-30/-60/30)
    0.012356   % (15/-15/15/-15)
    0.0072     % (30/-30/30/-30)
    0.0020     % (45/-45/45/-45)
    0.0056     % (60/-60/60/-60)
    0.011518    % (75/-75/75/-75)
];

eps_norm = eps / strain_limit;

% ============================================================
% 1) FIRST PLOT: only the 5 angle-ply cases (90/±theta/90)
% ============================================================
idx1 = 1:6;

Y1 = [TW(idx1), PK(idx1), eps_norm(idx1)];
labels1 = layups(idx1);

figure('Color','w','Position',[200 120 1050 520]);
b1 = bar(Y1, 'grouped'); grid on; box on;
ylim([0, max([Y1(:); 1]) * 1.15]);

set(gca, 'XTick', 1:numel(idx1), 'XTickLabel', labels1, 'FontSize', 11);
xtickangle(20);

ylabel('Failure index / Normalized strain');
title(sprintf('Multidirectional Laminate for Sphere'));

legend({'Tsai-Wu','Puck','Max Strain index'}, 'Location','northeast');

% failure/limit line at 1
hold on;
yline(1, '--', 'LineWidth', 1.2, HandleVisibility='off');
hold off;

% optional: value labels on bars
for k = 1:numel(b1)
    xt = b1(k).XEndPoints;
    yt = b1(k).YEndPoints;
    vals = b1(k).YData;
    text(xt, yt, compose('%.3f', vals), 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize', 9);
end

% ============================================================
% 2) SECOND PLOT: remaining 6 layups
% ============================================================
idx2 = 7:11;

Y2 = [TW(idx2), PK(idx2), eps_norm(idx2)];
labels2 = layups(idx2);

figure('Color','w','Position',[200 120 1050 520]);
b2 = bar(Y2, 'grouped'); grid on; box on;
ylim([0, max([Y2(:); 1]) * 1.15]);

set(gca, 'XTick', 1:numel(idx2), 'XTickLabel', labels2, 'FontSize', 11);
xtickangle(20);

ylabel('Failure index / Normalized strain');
title(sprintf('Angle ply Laminate for Sphere'));

legend({'Tsai-Wu','Puck','Max Strain index'}, 'Location','northeast');

hold on;
yline(1, '--', 'LineWidth', 1.2, 'HandleVisibility','off');
hold off;

% optional: value labels on bars
for k = 1:numel(b2)
    xt = b2(k).XEndPoints;
    yt = b2(k).YEndPoints;
    vals = b2(k).YData;
    text(xt, yt, compose('%.3f', vals), 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize', 9);
end

%% ============================================================
% n = 0.45
% ============================================================
strain_limit = 0.005;

layups = {
    '(90/15/-15/90)'
    '(90/30/-30/90)'
    '(90/45/-45/90)'
    '(90/60/-60/90)'
    '(90/75/-75/90)'
    '(60/-30/-60/30)'
    '(15/-15/15/-15)'
    '(30/-30/30/-30)'
    '(45/-45/45/-45)'
    '(60/-60/60/-60)'
    '(75/-75/75/-75)'
};

TW = [
    0.2266    % (90/15/-15/90)
    0.2436    % (90/30/-30/90)
    0.3009    % (90/45/-45/90)
    0.4578    % (90/60/-60/90)
    0.6335    % (90/75/-75/90)
    0.2771    % (60/-30/-60/30)
    1.0809    % (15/-15/15/-15)
    0.7064    % (30/-30/30/-30)
    0.3044    % (45/-45/45/-45)
    0.3272    % (60/-60/60/-60)
    0.5802    % (75/-75/75/-75)
];

PK = [
    0.1504    % (90/15/-15/90)
    0.1623    % (90/30/-30/90)
    0.2578    % (90/45/-45/90)
    0.3759    % (90/60/-60/90)
    0.5441    % (90/75/-75/90)
    0.2099    % (60/-30/-60/30)
    1.0510    % (15/-15/15/-15)
    0.6139    % (30/-30/30/-30)
    0.2100    % (45/-45/45/-45)
    0.2868    % (60/-60/60/-60)
    0.4838    % (75/-75/75/-75)
];

eps = [
    0.000726   % (90/15/-15/90)
    0.000989   % (90/30/-30/90)
    0.001300   % (90/45/-45/90)
    0.001883   % (90/60/-60/90)
    0.002343   % (90/75/-75/90)
    0.000917   % (60/-30/-60/30)
    0.004400   % (15/-15/15/-15)
    0.002500   % (30/-30/30/-30)
    0.000790   % (45/-45/45/-45)
    0.001288   % (60/-60/60/-60)
    0.002050   % (75/-75/75/-75)
];

eps_norm = eps / strain_limit;

% ============================================================
% 1) FIRST PLOT: only the 5 angle-ply cases (90/±theta/90)
% ============================================================
idx1 = 1:6;

Y1 = [TW(idx1), PK(idx1), eps_norm(idx1)];
labels1 = layups(idx1);

figure('Color','w','Position',[200 120 1050 520]);
b1 = bar(Y1, 'grouped'); grid on; box on;
ylim([0, max([Y1(:); 1]) * 1.15]);

set(gca, 'XTick', 1:numel(idx1), 'XTickLabel', labels1, 'FontSize', 11);
xtickangle(20);

ylabel('Failure index / Normalized strain');
title(sprintf('Multidirectional Laminate for Super Ellipsoid with n_1 = 0.45'));

legend({'Tsai-Wu','Puck','Max Strain index'}, 'Location','northeast');

% failure/limit line at 1
hold on;
yline(1, '--', 'LineWidth', 1.2, HandleVisibility='off');
hold off;

% optional: value labels on bars
for k = 1:numel(b1)
    xt = b1(k).XEndPoints;
    yt = b1(k).YEndPoints;
    vals = b1(k).YData;
    text(xt, yt, compose('%.3f', vals), 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize', 9);
end

% ============================================================
% 2) SECOND PLOT: remaining 6 layups
% ============================================================
idx2 = 7:11;

Y2 = [TW(idx2), PK(idx2), eps_norm(idx2)];
labels2 = layups(idx2);

figure('Color','w','Position',[200 120 1050 520]);
b2 = bar(Y2, 'grouped'); grid on; box on;
ylim([0, max([Y2(:); 1]) * 1.15]);

set(gca, 'XTick', 1:numel(idx2), 'XTickLabel', labels2, 'FontSize', 11);
xtickangle(20);

ylabel('Failure index / Normalized strain');
title(sprintf('Angle ply Laminate for Super Ellipsoid with n_1 = 0.45'));

legend({'Tsai-Wu','Puck','Max Strain index'}, 'Location','northeast');

hold on;
yline(1, '--', 'LineWidth', 1.2, 'HandleVisibility','off');
hold off;

% optional: value labels on bars
for k = 1:numel(b2)
    xt = b2(k).XEndPoints;
    yt = b2(k).YEndPoints;
    vals = b2(k).YData;
    text(xt, yt, compose('%.3f', vals), 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize', 9);
end
T  = [293 193  93];   % K
Xt = [3298 3183 3168]; % MPa
Xc = [1753 1732 1715]; % MPa
Yt = [  50   50   53]; % MPa
Yc = [  96  148  259]; % MPa
S12= [  81  126  169]; % MPa

Tq = 20; % target temperature (K)

% Use shape-preserving cubic interpolation
Xt_20  = pchip(T, Xt,  Tq);
Xc_20  = pchip(T, Xc,  Tq);
Yt_20  = pchip(T, Yt,  Tq);
Yc_20  = pchip(T, Yc,  Tq);
S12_20 = pchip(T, S12, Tq);

fprintf('Interpolated (pchip extrapolated) strengths at T = %g K:\n', Tq);
fprintf('Xt  = %.2f MPa\n', Xt_20);
fprintf('Xc  = %.2f MPa\n', Xc_20);
fprintf('Yt  = %.2f MPa\n', Yt_20);
fprintf('Yc  = %.2f MPa\n', Yc_20);
fprintf('S12 = %.2f MPa\n', S12_20);

Results = table(Tq, Xt_20, Xc_20, Yt_20, Yc_20, S12_20, ...
    'VariableNames', {'T_K','Xt_MPa','Xc_MPa','Yt_MPa','Yc_MPa','S12_MPa'});
disp(Results)

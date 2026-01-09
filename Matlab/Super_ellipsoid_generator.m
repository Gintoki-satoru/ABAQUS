% Super-ellipsoid parameters
a = 100;
b = 100;
c = 100;
thick = 1;
epsilon1 = 1;
epsilon2 = 1;

signed_power = @(base, exp) sign(base) .* (abs(base).^exp);

% Parameter grid
[u, v] = meshgrid(linspace(-pi, pi, 80), linspace(-pi/2, pi/2, 40));

% Super-ellipsoid equations
x = a * signed_power(cos(v), epsilon2) .* signed_power(cos(u), epsilon1);
y = b * signed_power(cos(v), epsilon2) .* signed_power(sin(u), epsilon1);
z = c * signed_power(sin(v), epsilon2);

% Plot
figure;
surf(z, x, y, 'FaceColor', 'texturemap', 'EdgeColor', 'texturemap');
axis equal;
xlabel('Z'); ylabel('X'); zlabel('Y');
title('Super-Ellipsoid');
grid on;
camlight; lighting gouraud;

%%
figure;
surf(z, x, y, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
axis equal;
axis off;              % ← hides axes, ticks, and labels
camlight; lighting gouraud;
%% Parameters for animation
a = 100;
b = 100;
c = 150;

% Range of exponents to animate
n_forward  = linspace(2, 0.2, 120);
n_backward = linspace(0.2, 2, 120);
n_vals = [n_forward n_backward(2:end)];

% Signed power helper
signed_power = @(base, exp) sign(base) .* (abs(base).^exp);

% Parameter grid (resolution controls smoothness)
[u, v] = meshgrid(linspace(-pi, pi, 120), linspace(-pi/2, pi/2, 60));

% Create video writer
vfile = VideoWriter('supern!.mp4','MPEG-4');
vfile.FrameRate = 30;
vfile.Quality   = 95;
open(vfile);

% Figure setup
fig = figure('Color', 'white');
axis equal; axis off;
view(75, 75);

% Animation loop
for k = 1:length(n_vals)

    n1 = n_vals(k);
    n2 = 1;

    % Parametric equations
    x = a * signed_power(cos(v), n1) .* signed_power(cos(u), n2);
    y = b * signed_power(cos(v), n1) .* signed_power(sin(u), n2);
    z = c * signed_power(sin(v), n1);

    % Plot
    % fig = figure('Color', '#004e73');
    % set(gca, 'Color', '#004e73');
    surf(z, x, y, 'FaceColor','interp','EdgeColor','none');
    axis equal; axis off;
    camlight headlight; lighting gouraud;

    % title(sprintf('Super-Ellipsoid Morphing (n1 = n2 = %.2f)', n1), ...
    %     'FontSize', 14);

    % Optional slow rotation
    % view(45 + 0.5*k, 25);

    % Capture frame
    frame = getframe(fig);
    writeVideo(vfile, frame);
end

close(vfile);

disp('Animation complete: SuperEllipsoid_n1n2_morph.mp4');

%% Sphere in a bounding box and save as image

% --- Parameters ---
r = 100;                        % sphere radius
center = [0 0 0];             % sphere center
box_scale = 1;              % box a bit larger than sphere

% --- Create sphere geometry ---
[XS, YS, ZS] = sphere(100);   % unit sphere
XS = r * XS + center(1);
YS = r * YS + center(2);
ZS = r * ZS + center(3);

% --- Bounding box (cube) coordinates ---
L = box_scale * r;
c = center;

corners = [ ...
    c + [-L -L -L];  % 1
    c + [ L -L -L];  % 2
    c + [ L  L -L];  % 3
    c + [-L  L -L];  % 4
    c + [-L -L  L];  % 5
    c + [ L -L  L];  % 6
    c + [ L  L  L];  % 7
    c + [-L  L  L]]; % 8

edges = [ ...
    1 2; 2 3; 3 4; 4 1; ...  % bottom square
    5 6; 6 7; 7 8; 8 5; ...  % top square
    1 5; 2 6; 3 7; 4 8];     % vertical edges

% --- Plot ---
fig = figure('Color','white');
ax = axes('Parent', fig);
hold(ax, 'on');

% Sphere
surf(ax, XS, YS, ZS, 'FaceColor',"#FDCA00", 'EdgeColor','none');
camlight; lighting gouraud;

% Bounding box
for k = 1:size(edges,1)
    p1 = corners(edges(k,1),:);
    p2 = corners(edges(k,2),:);
    plot3(ax, [p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], ...
        'k-', 'LineWidth', 1.5);
end

axis equal off
view(35, 25);

% --- Save as image ---
% exportgraphics(fig, 'sphere_bbox.png', 'Resolution', 300);
% disp('Saved image as sphere_bbox.png');

%% === Geometry Parameters ===
r = 1;                 % radius of cylinder + spherical heads
Lc = 2;                % cylindrical section length
% Total vessel length = cylinder length + 2 radii (spherical heads)
L_total = Lc + 2*r;

% === Parameter Grid ===
N = 80;
theta = linspace(0, 2*pi, N)';
z_cyl = linspace(-Lc/2, Lc/2, N);

% Cylinder coordinates
Xc = r * cos(theta) * ones(1, N);
Yc = r * sin(theta) * ones(1, N);
Zc = ones(N,1) * z_cyl;

% Spherical end caps (left & right)
[phi, th] = meshgrid(linspace(0, pi, N), linspace(0, 2*pi, N));

% Left spherical head (center at -Lc/2)
Xs1 = r * sin(phi) .* cos(th);
Ys1 = r * sin(phi) .* sin(th);
Zs1 = -Lc/2 - r * cos(phi);

% Right spherical head (center at +Lc/2)
Xs2 = r * sin(phi) .* cos(th);
Ys2 = r * sin(phi) .* sin(th);
Zs2 = +Lc/2 + r * cos(phi);
% === Bounding Box Dimensions ===
box_min = [-r, -r, -L_total/2];
box_max = [+r, +r, +L_total/2];

% Extract corner points
xmn = box_min(1); ymn = box_min(2); zmn = box_min(3);
xmx = box_max(1); ymx = box_max(2); zmx = box_max(3);

corners = [
    xmn ymn zmn;  xmx ymn zmn;  xmx ymx zmn;  xmn ymx zmn;  % bottom
    xmn ymn zmx;  xmx ymn zmx;  xmx ymx zmx;  xmn ymx zmx   % top
];

edges = [
    1 2; 2 3; 3 4; 4 1; ...
    5 6; 6 7; 7 8; 8 5; ...
    1 5; 2 6; 3 7; 4 8
];

% === Plotting ===
fig = figure('Color', 'white');
hold on;

% Cylinder
surf(Xc, Yc, Zc, 'FaceColor',"#AFCC50", 'EdgeColor','none');

% Spherical heads
surf(Xs1, Ys1, Zs1, 'FaceColor',"#AFCC50", 'EdgeColor','none');
surf(Xs2, Ys2, Zs2, 'FaceColor',"#AFCC50", 'EdgeColor','none');

% Bounding box
for k = 1:size(edges,1)
    p1 = corners(edges(k,1), :);
    p2 = corners(edges(k,2), :);
    plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], ...
        'k-', 'LineWidth', 1.3);
end

axis equal off
view(35, 20);
camlight; lighting gouraud;

%% ===== Super-Ellipsoid Parameters =====
a = 1;
b = 1;
c = 1;

n1 = 0.1;   % exponent along u
n2 = 1;   % exponent along v

signed_power = @(base, exp) sign(base) .* (abs(base).^exp);

% ===== Parameter Grid =====
N = 120;
[u, v] = meshgrid(linspace(-pi, pi, N), linspace(-pi/2, pi/2, N));

% ===== Super-Ellipsoid Parametric Surface =====
X = a * signed_power(cos(v), n1) .* signed_power(cos(u), n2);
Y = b * signed_power(cos(v), n1) .* signed_power(sin(u), n2);
Z = c * signed_power(sin(v), n1);

% ===== Compute Bounding Box (tight) =====
xmin = min(X(:));  xmax = max(X(:));
ymin = min(Y(:));  ymax = max(Y(:));
zmin = min(Z(:));  zmax = max(Z(:));

corners = [
    xmin ymin zmin;
    xmax ymin zmin;
    xmax ymax zmin;
    xmin ymax zmin;
    xmin ymin zmax;
    xmax ymin zmax;
    xmax ymax zmax;
    xmin ymax zmax
];

edges = [
    1 2; 2 3; 3 4; 4 1;
    5 6; 6 7; 7 8; 8 5;
    1 5; 2 6; 3 7; 4 8
];

% ===== Plotting =====
fig = figure('Color','white');
hold on;

% Super-ellipsoid surface
surf(X, Y, Z, 'FaceColor',"#EE7A34", 'EdgeColor','none');
camlight headlight; lighting gouraud;

% Bounding box
for k = 1:size(edges,1)
    p1 = corners(edges(k,1),:);
    p2 = corners(edges(k,2),:);
    plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], ...
          'k-', 'LineWidth', 1.4);
end

axis equal off
view(35, 20);


%%

a = 1;
b = 1;
n_list = [0.5 1 2 3 4];

% helper: signed power
spow = @(x,p) sign(x) .* (abs(x).^p);

t = linspace(0, 2*pi, 1000);

figure; hold on; box on;

for n = n_list
    % Lamé curve: (|x/a|)^n + (|y/b|)^n = 1
    % Parametric form:
    x = a * spow(cos(t), 2/n);
    y = b * spow(sin(t), 2/n);

    plot(x, y, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('n = %.1f', n));
end

axis equal;
xlabel('a');
ylabel('b');
title('Lamé curves for a = b = 1 and different n');
legend('Location','bestoutside');
grid on;

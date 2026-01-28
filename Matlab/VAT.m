function bendfree_superellipsoid_demo()
% Bend-free VAT design for super-ellipsoid of revolution (stress-focused)
% Outputs:
%   - Delta(phi): curvature ratio r_theta/r_phi
%   - (V1*, V2*): tailored lamination parameters along phi
%   - theta(phi): equivalent tow angle distribution (deg)
%
% Reference: Daghighi et al., "Bend-free design of super ellipsoids of revolution
% composite pressure vessels", Composite Structures 245 (2020) 112283.
%
% Notes:
% - This script focuses on stress/bending suppression (bend-free), not failure.
% - Uses numerical curvature from meridian curve (robust for any N, K).

%% ---------------- USER INPUTS ----------------
a = 114.24;          % semi-axis in radial direction (mm)
K = 7;         % aspect ratio b/a
b = K*a;          % semi-axis in axial direction (mm)
N = 2/0.7;          % shape parameter (Lamé exponent)
nPhi = 181;       % number of phi samples (0..90 deg)

% UD lamina properties in MPa - 20K car epx
E11  = 173930;    % [MPa]
E22  = 15690;     % [MPa]
G12 = 11244;      % [MPa]
nu12 = 0.433;

% Choose feasible-region "lower boundary" model.
% For balanced/symmetric feasible region (typical 2D map), a common lower
% boundary shape is: V2 = -sqrt((V1^2 + 1)/2). This matches the "parabola"
% style boundary seen in the paper's feasible region figure.
useLowerBoundaryParabola = true;

% ---------------- COMPUTE MATERIAL INVARIANTS Ui ----------------
U = lamina_invariants_U(E11, E22, G12, nu12); % returns struct U1..U4

% ---------------- SAMPLE PHI AND GEOMETRY ----------------
phi = linspace(0, 90, nPhi);       % degrees
phiRad = deg2rad(phi);

% Meridian curve of superellipse in r-z plane:
% (r/a)^N + (z/b)^N = 1 for first quadrant.
% Parametrization:
%   r = a * cos(phi)^(2/N)
%   z = b * sin(phi)^(2/N)
% This is a standard superellipse parametrization for quadrant.
[r, z] = superellipse_meridian(a, b, N, phiRad);

% Curvature ratio Delta(phi) = r_theta / r_phi computed from principal curvatures
Delta = curvature_ratio_delta(r, z);

% ---------------- BEND-FREE LINE IN V1-V2 SPACE (Eq. 14 -> linear) ----------------
% From the parsed equation in the paper, Eq. (14) can be rearranged to:
%   A1*V1 + A2*V2 + A0 = 0  (for each phi, i.e., each Delta)
% which is a straight line in (V1, V2).
[A0, A1, A2] = bendfree_line_coeffs(Delta, U);

% Define line: V2 = m*V1 + c
m = -A1 ./ A2;
c = -A0 ./ A2;

% ---------------- INTERSECT WITH LOWER FEASIBLE BOUNDARY ----------------
V1star = nan(size(phi));
V2star = nan(size(phi));

for i = 1:numel(phi)
    if ~isfinite(m(i)) || ~isfinite(c(i))
        continue;
    end

    % Solve intersection in V1 within [-1, 1]
    % g(V1) = V2_line(V1) - V2_lower(V1) = 0
    g = @(V1) (m(i)*V1 + c(i)) - V2_lower_boundary(V1, useLowerBoundaryParabola);

    % Try to find a sign change bracket
    V1grid = linspace(-1, 1, 2001);
    gvals = arrayfun(g, V1grid);
    idx = find(isfinite(gvals), 1, 'first');
    if isempty(idx), continue; end

    % Find intervals where g changes sign
    s = sign(gvals);
    s(~isfinite(s)) = 0;
    changeIdx = find(s(1:end-1).*s(2:end) <= 0 & s(1:end-1)~=0 & s(2:end)~=0, 1, 'first');

    if isempty(changeIdx)
        % No intersection found in [-1,1] with this boundary model
        continue;
    end

    xL = V1grid(changeIdx);
    xU = V1grid(changeIdx+1);

    % Root solve
    V1sol = fzero(g, [xL, xU]);
    V2sol = m(i)*V1sol + c(i);

    % Clip to physical bounds
    if abs(V1sol) <= 1+1e-6 && abs(V2sol) <= 1+1e-6
        V1star(i) = max(-1, min(1, V1sol));
        V2star(i) = max(-1, min(1, V2sol));
    end
end

% ---------------- MAP (V1*, V2*) -> EQUIVALENT TOW ANGLE theta(phi) ----------------
% Use a single equivalent angle theta in [0,90] such that:
%   V1(theta)=cos(2 theta), V2(theta)=cos(4 theta)
% Then minimize (V1-V1*)^2 + (V2-V2*)^2.
thetaDeg = nan(size(phi));
for i = 1:numel(phi)
    if ~isfinite(V1star(i)) || ~isfinite(V2star(i)), continue; end
    obj = @(th) (cosd(2*th)-V1star(i)).^2 + (cosd(4*th)-V2star(i)).^2;
    thetaDeg(i) = fminbnd(obj, 0, 90);
end

% ---------------- PLOTS (useful for stress comparison work) ----------------
figure; plot(phi, Delta, 'LineWidth', 1.5);
xlabel('\phi (deg)'); ylabel('\Delta(\phi)=r_\theta/r_\phi');
title(sprintf('Curvature ratio \\Delta along super-ellipsoid (N=%.2f, K=%.2f)', N, K));
grid on;

figure; plot(phi, V1star, 'LineWidth', 1.5); hold on;
plot(phi, V2star, 'LineWidth', 1.5);
xlabel('\phi (deg)'); ylabel('Tailored lamination parameters');
legend('V_1^*','V_2^*','Location','best');
title('Bend-free tailored lamination parameters along \phi');
grid on;

figure; plot(phi, thetaDeg, 'LineWidth', 1.5);
xlabel('\phi (deg)'); ylabel('\theta(\phi) (deg)');
title('Equivalent VAT tow angle distribution for bend-free state');
grid on;

%% ---------------- EXPORT TABLE ----------------
T = table(phi(:), Delta(:), V1star(:), V2star(:), thetaDeg(:), ...
    'VariableNames', {'phi_deg','Delta','V1_star','V2_star','theta_deg'});
disp(T(1:10,:));
assignin('base','bendfree_table',T);

end

%% ========================================================================
function [r, z] = superellipse_meridian(a, b, N, phiRad)
% Meridian superellipse (first quadrant) parametrization:
% r = a*cos(phi)^(2/N), z = b*sin(phi)^(2/N)
c = cos(phiRad);
s = sin(phiRad);

% Ensure non-negative (first quadrant)
c = max(c, 0);
s = max(s, 0);

r = a * (c).^(2./N);
z = b * (s).^(2./N);

% Avoid r=0 exactly at pole for curvature calc stability (tiny epsilon)
r(1) = max(r(1), 1e-9*a);
end

function Delta = curvature_ratio_delta(r, z)
% Compute principal curvature radii for surface of revolution from meridian curve.
% Meridian curve: (r(φ), z(φ))
% Let s = arc length parameter implicitly via derivatives wrt φ index spacing.
% Principal curvatures:
%   k_mer = (r' z'' - z' r'') / (r'^2 + z'^2)^(3/2)
%   k_circ = z' / (r * sqrt(r'^2 + z'^2))
% Radii:
%   r_phi = 1/|k_mer|, r_theta = 1/|k_circ|
% Delta = r_theta/r_phi = |k_mer|/|k_circ|

% Use numerical derivatives wrt index; scale cancels in ratio.
dr = gradient(r);
dz = gradient(z);
d2r = gradient(dr);
d2z = gradient(dz);

den = (dr.^2 + dz.^2).^(3/2);
k_mer = (dr.*d2z - dz.*d2r) ./ max(den, eps);

k_circ = dz ./ max(r .* sqrt(dr.^2 + dz.^2), eps);

Delta = abs(k_mer) ./ max(abs(k_circ), eps);

% Clean near singular points
Delta(~isfinite(Delta)) = nan;
end

function U = lamina_invariants_U(E11, E22, G12, nu12)
% Compute reduced stiffness Q and then invariants U1..U4 as used in bend-free formulation.
% For plane stress orthotropic lamina:
nu21 = nu12 * (E22/E11);
den = 1 - nu12*nu21;

Q11 = E11/den;
Q22 = E22/den;
Q12 = nu12*E22/den;
Q66 = G12;

% Common Tsai/Hahn invariants used in lamination parameter forms (A-matrix).
% One consistent set (widely used with lamination parameters):
U1 = (3*Q11 + 3*Q22 + 2*Q12 + 4*Q66)/8;
U2 = (Q11 - Q22)/2;
U3 = (Q11 + Q22 - 2*Q12 - 4*Q66)/8;
U4 = (Q11 + Q22 + 6*Q12 - 4*Q66)/8;

U.U1 = U1; U.U2 = U2; U.U3 = U3; U.U4 = U4;
end

function [A0, A1, A2] = bendfree_line_coeffs(Delta, U)
% Rearranged bend-free governing equation (Eq. 14 in the paper) into:
%   A1*V1 + A2*V2 + A0 = 0
%
% Based on the expression shown in the parsed PDF:
%   -Δ^2 + 2Δ*(U1 - U4 + U2 V1 + U3 V2)/(U1 + U2 V1 + U3 V2)
%         + (2U4 - U1 - 3U3 V2 + U2 V1)/(U1 + U2 V1 + U3 V2) = 0
%
% Multiply by D = U1 + U2 V1 + U3 V2 and collect terms.

U1 = U.U1; U2 = U.U2; U3 = U.U3; U4 = U.U4;

A1 = (-Delta.^2).*U2 + 2.*Delta.*U2 + U2;
A2 = (-Delta.^2).*U3 + 2.*Delta.*U3 - 3.*U3;
A0 = (-Delta.^2).*U1 + 2.*Delta.*(U1 - U4) + (2.*U4 - U1);
end

function V2low = V2_lower_boundary(V1, useParabola)
% Lower feasible boundary model (balanced/symmetric lamination parameter space)
% If useParabola = true:
%   V2 = -sqrt((V1^2 + 1)/2)   (parabola-like lower boundary)
% Else:
%   fallback: V2 = -1 (crude bound)

if useParabola
    V2low = -sqrt((V1.^2 + 1)./2);
else
    V2low = -1 + 0*V1;
end
end

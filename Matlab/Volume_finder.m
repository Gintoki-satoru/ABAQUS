% --- Target volumes from the reference tank ---
V_inner_target = 0.0587;    % m³
V_wall_target  = 0.0017496;    % m³
V_outer_target = V_inner_target + V_wall_target;

% --- Superellipsoid exponents ---
n1 = 0.6;
n2 = 1.4;

% --- Aspect ratio (a:b:c) ---
ratio = [1, 1, 3];

% --- Initial guess: [a, t] ---
x0 = [0.3, 0.004];  % a ~ 0.3 m, t ~ 4 mm

% --- Solve for a and t ---
opts = optimoptions('fsolve','Display','iter','TolFun',1e-9);
x_sol = fsolve(@(x) volume_error_ratio(x, ratio, V_inner_target, V_outer_target, n1, n2), x0, opts);

a = x_sol(1);
t = x_sol(2);
b = (ratio(2)/ratio(1)) * a;
c = (ratio(3)/ratio(1)) * a;

fprintf('\nEquivalent superellipsoid:\n');
fprintf(' a = %.4f m\n b = %.4f m\n c = %.4f m\n thickness = %.5f m\n', a, b, c, t);

% ================================================================
% --- Helper functions -------------------------------------------
% ================================================================

function F = volume_error_ratio(x, ratio, V_inner_target, V_outer_target, n1, n2)
    a = x(1); 
    t = x(2);
    b = (ratio(2)/ratio(1))*a;
    c = (ratio(3)/ratio(1))*a;

    % Outer dimensions preserve same ratio
    a_o = a + t;
    b_o = b + t;
    c_o = c + t;

    % Compute inner & outer volumes
    V_inner = superellipsoid_volume(a, b, c, n1, n2);
    V_outer = superellipsoid_volume(a_o, b_o, c_o, n1, n2);

    % Two equations: match inner and wall volumes
    F = [V_inner - V_inner_target;
         (V_outer - V_inner) - (V_outer_target - V_inner_target)];
end

function V = superellipsoid_volume(a, b, c, n1, n2)
    e1 = 2/n1;
    e2 = 2/n2;
    V = 8 * a * b * c * (gamma(1 + 1/e1)^2 * gamma(1 + 1/e2)) / ...
        (gamma(1 + 2/e1) * gamma(1 + (1/e2 + 2/e1)));
end

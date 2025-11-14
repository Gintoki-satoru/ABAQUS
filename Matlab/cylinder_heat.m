r_inner = 0.85;      % m
t_liner = 0.00326;   % m
t_ins   = 0.016;     % m
t_outer = 0.005;     % m
L       = 1.95;      % m
h_out   = 10;        % W/m2K
k_liner = 30.6;      % W/mK
k_ins   = 3.03e-5;    % W/mK
T_inside  = 20;      % K
T_ambient = 300;     % K

Q = total_heat_flow(r_inner, t_liner, t_ins, t_outer, L, h_out, k_liner, k_ins, T_inside, T_ambient)
r_liner_out = r_inner + t_liner;
A_dome_liner = 4*pi*r_liner_out^2;
A_cyl_liner = 2*pi*r_liner_out*L;
A_liner = (A_dome_liner + A_cyl_liner)*10^6;
heq = Q/(A_liner*(T_ambient - T_inside))

function Q = total_heat_flow(r_inner, t_liner, t_ins, t_outer, L, h_out, k_liner, k_ins, T_inside, T_ambient)
% TOTAL_HEAT_FLOW  Compute total heat leak Q through cylindrical tank
% including liner, insulation, outer wall, and convection.
%
% Inputs:
%   r_inner     - inner radius of liner (m)
%   t_liner     - liner thickness (m)
%   t_ins       - insulation thickness (m)
%   t_outer     - outer tank thickness (m)
%   L           - length of cylindrical section (m)
%   h_out       - convective heat transfer coeff outside (W/m^2K)
%   k_liner     - thermal conductivity of liner (W/mK)
%   k_ins       - thermal conductivity of insulation (W/mK)
%   T_inside    - inner temperature (K)
%   T_ambient   - ambient temperature (K)
%
% Output:
%   Q           - total heat flow (W)

% ---- Compute Radii ----
r_liner_out = r_inner + t_liner;
r_ins_out   = r_liner_out + t_ins;
r_outer     = r_ins_out + t_outer;

% ---- Cylindrical Surfaces ----
A_cyl_out   = 2*pi*r_outer*L;

% ---- Dome Areas (two hemispheres = one sphere) ----
A_dome_out   = 4*pi*r_outer^2;

% ---- Liner Resistance (cyl + domes in parallel) ----
R_liner_cyl = log(r_liner_out / r_inner) / (2*pi*k_liner*L);
R_liner_dome = (1/(4*pi*k_liner)) * (1/r_inner - 1/r_liner_out);

R_liner = (R_liner_cyl * R_liner_dome) / (R_liner_cyl + R_liner_dome);

% ---- Insulation Resistance (cyl + domes in parallel) ----
R_ins_cyl = log(r_ins_out / r_liner_out) / (2*pi*k_ins*L);
R_ins_dome = (1/(4*pi*k_ins)) * (1/r_liner_out - 1/r_ins_out);

R_ins = (R_ins_cyl * R_ins_dome) / (R_ins_cyl + R_ins_dome);

% ---- Outer convection ----
A_outer_total = A_cyl_out + A_dome_out;
R_conv = 1 / (h_out * A_outer_total);

% ---- Total resistance ----
R_total = R_liner + R_ins + R_conv;

% ---- Heat flow ----
Q = (T_ambient - T_inside) / R_total;
end




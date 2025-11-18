% ---- INPUTS ----
r_inner   = 0.305/2;       % m
t_liner   = 0.002;    % m
t_ins     = 0.016;      % m
t_outer   = 0.002;      % m
L         = 0.6;       % m

h_out     = 10;         % W/m^2K
k_liner   = 30.6;       % W/mK
k_ins     = 3.03e-5;    % W/mK  % 0.0303 mW/mK
k_outer   = 30.6;       % W/mK

T_inside  = 20;         % K
T_ambient = 300;        % K

Q = total_heat_flow(r_inner, t_liner, t_ins, t_outer, L, ...
                    h_out, k_liner, k_ins, k_outer, ...
                    T_inside, T_ambient)

% Equivalent heat transfer coefficient based on LINER outer area:
r_liner_out = r_inner + t_liner;
A_dome_liner = 4*pi*r_liner_out^2;
A_cyl_liner  = 2*pi*r_liner_out*L;
A_liner      = (A_dome_liner + A_cyl_liner) * 1e6;  % [mm^2]
heq          = Q / (A_liner * (T_ambient - T_inside))  % [W/mm^2K]


function Q = total_heat_flow(r_inner, t_liner, t_ins, t_outer, L, ...
                             h_out, k_liner, k_ins, k_outer, ...
                             T_inside, T_ambient)
% TOTAL_HEAT_FLOW  Compute total heat leak Q through cylindrical tank
% including liner, insulation, OUTER WALL, and convection.
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
%   k_outer     - thermal conductivity of outer tank wall (W/mK)
%   T_inside    - inner temperature (K)
%   T_ambient   - ambient temperature (K)
%
% Output:
%   Q           - total heat flow (W)

% ---- Radii ----
r_liner_out = r_inner + t_liner;
r_ins_out   = r_liner_out + t_ins;
r_outer     = r_ins_out + t_outer;

% ---- Outer surfaces (for convection) ----
A_cyl_out  = 2*pi*r_outer*L;
A_dome_out = 4*pi*r_outer^2;
A_outer_total = A_cyl_out + A_dome_out;

% ---- Liner resistance (cylinder + domes in parallel) ----
R_liner_cyl  = log(r_liner_out / r_inner) / (2*pi*k_liner*L);
R_liner_dome = (1/(4*pi*k_liner)) * (1/r_inner - 1/r_liner_out);
R_liner      = (R_liner_cyl * R_liner_dome) / (R_liner_cyl + R_liner_dome);

% ---- Insulation resistance (cylinder + domes in parallel) ----
R_ins_cyl  = log(r_ins_out / r_liner_out) / (2*pi*k_ins*L);
R_ins_dome = (1/(4*pi*k_ins)) * (1/r_liner_out - 1/r_ins_out);
R_ins      = (R_ins_cyl * R_ins_dome) / (R_ins_cyl + R_ins_dome);

% ---- Outer wall resistance (cylinder + domes in parallel) ----
R_outer_cyl  = log(r_outer / r_ins_out) / (2*pi*k_outer*L);
R_outer_dome = (1/(4*pi*k_outer)) * (1/r_ins_out - 1/r_outer);
R_outer      = (R_outer_cyl * R_outer_dome) / (R_outer_cyl + R_outer_dome);

% ---- Outer convection ----
R_conv = 1 / (h_out * A_outer_total);

% ---- Total resistance (series) ----
R_total = R_liner + R_ins + R_outer + R_conv;

% ---- Heat flow ----
Q = 1.13*1.3*(T_ambient - T_inside) / R_total;
end

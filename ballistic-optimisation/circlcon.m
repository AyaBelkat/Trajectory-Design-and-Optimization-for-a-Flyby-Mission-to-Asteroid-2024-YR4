function [c, ceq] = circlcon(x)
% =========================================================================
% Nonlinear constraints: ballistic Earth → asteroid transfer
%
% Purpose:
% Enforce basic feasibility for a purely ballistic Lambert arc from Earth
% at t0 to asteroid 2024 YR4 at tf. Constraints cap departure/arrival
% hyperbolic excess, impose a minimum TOF, and encourage a reasonable
% arrival approach geometry.
%
% Decision vector:
% x = [t0, tf] % ET seconds (TDB) past J2000
%
% Modeling notes:
% • Sun-centered two-body propagation (patched-conics at endpoints).
% • States from SPICE in J2000 with no aberration correction ('NONE').
% • Lambert arc choose the lowest-|v1| branch as a pragmatic energy proxy.
%
% Inequality constraints assembled in c (must satisfy c ≤ 0):
% 1) cos(60°) − cos(arrival_angle) ≤ 0 % arrival direction (sunward approach)
% 2) ||v∞dep|| − 3.873 ≤ 0 % cap departure v∞ (km/s)
% 3) ||v∞arr|| − 7 ≤ 0 % cap arrival v∞ (km/s)
% 4) (100 d) − TOF ≤ 0 % minimum time of flight
%
% Equality constraints: none (ceq = []).
% =========================================================================

t0 = x(1);
tf = x(2);

% --- Constants and GTO reference orbit (used only for consistency) --------
muS = 1.32712440018e11; % Sun gravitational parameter [km^3/s^2]
muE = 398600; % Earth gravitational parameter [km^3/s^2]
Re = 6378.137; % Earth equatorial radius [km]
rp = Re + 250; % GTO perigee radius [km]
ra = Re + 22500; % GTO apogee radius [km]
a = 0.5*(ra + rp); % GTO semi-major axis [km]

% --- States at endpoints (Sun/J2000 frame) --------------------------------
D_earth = cspice_spkezr('Earth', t0, 'J2000', 'NONE', 'Sun');
D_ast = cspice_spkezr('54509621', tf, 'J2000', 'NONE', 'Sun');

rE = D_earth(1:3);
rAst = D_ast(1:3);
vE = D_earth(4:6);
vAst = D_ast(4:6);

% --- Lambert solve Earth→asteroid -----------------------------------------
TOF = tf - t0;
[v1, v2] = lambert_solver(rE, rAst, TOF, muS);

% If multiple solutions exist, pick the branch with the smallest |v1|
if size(v1, 1) > 1
    v0 = norm(v1(1, :));
    v11 = norm(v1(2, :));
    v22 = norm(v1(3, :));
    if v0 < v11 && v0 < v22
        v1 = v1(1, :);
        v2 = v2(1, :);
    elseif v11 < v0 && v11 < v22
        v1 = v1(2, :);
        v2 = v2(2, :);
    else
        v1 = v1(3, :);
        v2 = v2(3, :);
    end
end

% --- Hyperbolic excess at departure and arrival ---------------------------
v_inf1 = v1(:) - vE; % departure v∞ (Earth-relative)
v_inf1_n = norm(v_inf1);
v_inf2 = v2(:) - vAst; % arrival v∞ (asteroid-relative)
v_inf2_n = norm(v_inf2);

% --- Arrival geometry  ------------------------------------------
dir_rAst = -rAst/norm(rAst); % anti-solar radial direction at arrival
dir_v_arr = -v_inf2/v_inf2_n; % opposite to arrival v∞ (approach direction)

% Calculate cosine of angle between approach direction and anti-solar direction
cos_arr_angle = dot(dir_v_arr, dir_rAst);

% Clamp to valid range for numerical stability
cos_arr_angle = max(-1, min(1, cos_arr_angle));

% --- Assemble inequality constraints c ≤ 0 --------------------------------
min_TOF = 100*86400;        % overall minimum TOF [s]
c = [];

% (1) Arrival approach geometry: prefer sunward approach (≤60° from −r̂)
% Constraint: cos(arrival_angle) ≥ cos(60°)
% Rearranged: cos(60°) - cos(arrival_angle) ≤ 0
c(end + 1) = cosd(60) - cos_arr_angle;

% (2) Cap on departure v∞ to reflect launcher/escape energy constraints:
c(end + 1) = v_inf1_n - 3.873;

% (3) Arrival v∞ cap: ||v∞arr|| ≤ 7 km/s
c(end + 1) = v_inf2_n - 7;

% (4) Minimum time of flight
c(end + 1) = (t0 - tf);

% --- No equality constraints ----------------------------------------------
ceq = [];

end

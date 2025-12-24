function J = cost_function_ballistic_velocity(x)

% =========================================================================
% COST_FUNCTION_BALLISTIC_VELOCITY
%
% Purpose:
%   Defines the scalar objective for a purely ballistic Earth→asteroid
%   transfer optimization. The cost reflects the launch injection Δv
%   (from a reference GTO) needed to achieve the required departure
%   hyperbolic excess velocity.
%
% Decision vector:
%   x = [t0, tf]
%       t0 : departure epoch (ET, TDB seconds past J2000)
%       tf : arrival   epoch (ET, TDB seconds past J2000)
%
% Modelling assumptions:
%   • Heliocentric two-body dynamics (Sun GM).
%   • Planetary/asteroid states from SPICE (Sun-centred J2000 frame).
%   • Patched-conic at Earth departure, using GTO injection cost proxy.
%   • No midcourse DSMs or gravity assists; purely ballistic.
%
% Objective:
%   J = Δv_inj
%   i.e. the injection Δv from GTO perigee.
%
% Debug output:
%   Prints epochs, v∞ at departure/arrival, and arrival phase angle.
% =========================================================================

t0 = x(1);
tf = x(2);

% --- Constants (Sun + Earth GM, Earth radius, GTO reference orbit) --------
muS = 1.32712440018e11;   % Sun gravitational parameter [km^3/s^2]
muE = 398600;             % Earth gravitational parameter [km^3/s^2]
Re  = 6378.137;           % Earth equatorial radius [km]
rp  = Re + 250;           % GTO perigee radius [km]
ra  = Re + 22500;         % GTO apogee radius [km]
a   = 0.5*(ra + rp);      % GTO semi-major axis [km]

% --- SPICE ephemerides for Earth (at t0) and asteroid (at tf) --------------
D_earth = cspice_spkezr('Earth', t0, 'J2000', 'NONE', 'Sun');
D_ast = cspice_spkezr('54509621',  tf, 'J2000', 'NONE', 'Sun');  % 2024 YR4

rE = D_earth(1:3);
rAst = D_ast(1:3);

vE = D_earth(4:6);
vAst = D_ast(4:6);

% --- Lambert arc: Earth→asteroid ------------------------------------------
TOF = tf - t0;
[v1, v2] = lambert_solver(rE, rAst, TOF, muS);

% If multiple Lambert branches exist, select the lowest |v1| branch
if size(v1,1) > 1
    v0 = norm(v1(1, :));
    v11 = norm(v1(2, :));
    v22 = norm(v1(3, :));
    if v0 < v11 && v0 < v22
    v1 = v1(1, :);
    v2 = v2(1, :);
    end
    if v11 < v0 && v11 < v22
    v1 = v1(2, :);
    v2 = v2(2, :);
    end
    if v22 < v0 && v22 < v11 
    v1 = v1(3, :);
    v2 = v2(3, :);
    end
end

% --- Hyperbolic excess velocities -----------------------------------------
v_inf1    = v1' - vE;               % departure v∞ wrt Earth [km/s]     
v_inf1_n  = norm(v_inf1);
v_inf2 = v2' - vAst;                % arrival v∞ wrt asteroid [km/s]
v_inf2_n = norm(v_inf2);

% --- Injection Δv from GTO perigee ----------------------------------------
% Δv0 = v_peri(hyper) − v_peri(ellipse)
delta_v_inj = sqrt(v_inf1_n^2 + 2*muE/rp) - sqrt(muE*(2/rp - 1/a));

% --- Arrival phase angle diagnostic (not used in J) -----------------------
dir_rAst = -rAst/norm(rAst);            % anti-radial direction at arrival
dir_v_arr = -v_inf2/v_inf2_n;           % opposite to arrival v∞
arr_phase = dot(dir_v_arr, dir_rAst);   % angle [rad]
arr_phase = acosd(arr_phase);           % angle [deg]

% --- Objective function ---------------------------------------------------
Total_delta_v = delta_v_inj;    % cost = injection Δv only

% --- Debug logging --------------------------------------------------------
debug_log = true;
if debug_log  % log while optimizing
    fprintf('t0=%s, tf=%s, deltav_inj=%.3f km/s, v_inf2=%.3f km/s\n, arr_phase=%.3f deg  \n', ...
        cspice_et2utc(t0,'C',3), ...
        cspice_et2utc(tf,'C',3), ...
        delta_v_inj, v_inf2_n, arr_phase);
end


J = Total_delta_v;


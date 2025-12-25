function [c, ceq] = circlcon_GA_dsm(x)

% =========================================================================
% Nonlinear constraints for an Earth GA+DSM heliocentric transfer
%
% Purpose (dissertation-style summary):
%   Enforces feasibility for a trajectory with one unpowered Earth gravity
%   assist (GA) and one deep-space manoeuvre (DSM) under patched-conics.
%   The decision vector is
%       x = [t0, TOF, eta1, eta2, theta_b, Rperi]
%   where:
%       t0      : departure epoch (ET seconds, TDB)
%       TOF     : total time of flight (s)
%       eta1    : fraction of TOF at which GA occurs (Earth→GA duration = eta1*TOF)
%       eta2    : fraction of post-GA leg when DSM occurs (GA→DSM = eta2*(1-eta1)*TOF)
%       theta_b : B-plane angle selector for the GA (rad)
%       Rperi   : Earth flyby periapsis radius (km) (>= safety margin)
%
%   Inequality constraints c <= 0 include:
%     • minimum segment durations (overall and per leg),
%     • maximum departure hyperbolic excess speed ||v∞,dep||,
%     • arrival-geometry pointing (v∞,arr approximately anti-solar).
%   No equality constraints are imposed
%   except for the Lambert/propagation relations handled inside the cost.
%
% Usage:
%   Called by fmincon via nonlcon= @circlcon_GA_dsm in the GA+DSM optimizer.
%
% Notes:
%   • States are queried from SPICE (cspice_spkezr) in the Sun-centered
%     inertial frame (J2000). Lambert arcs are solved in two-body Sun GM.
%   • The gravity assist is modeled via an instantaneous rotation of the
%     incoming v∞ at Earth pericenter (classical unpowered flyby model).
%   • DSM is a single impulsive Δv at a chosen fraction along the post-GA leg.
% =========================================================================


t0 = x(1);
TOF = x(2);
eta1 = x(3);
eta2 = x(4);
theta_b = x(5);
Rperi = x(6);

% -------------------- Constants and basic GTO figures ---------------------
muS = 1.32712440018e11;   % Sun's gravitational parameter [km^3/s^2]
muE = 398600;             % Earth's gravitational parameter [km^3/s^2]
Re  = 6378.137;           % Radius of Earth [km]

% Earth GTO is used to quantify the initial injection cost.
% (rp, ra) allow computing Δv at perigee to reach the required v∞ wrt Earth.
rp  = Re + 250;           % GTO perigee radius [km]
ra  = Re + 22500;         % GTO apogee radius [km]
a   = 0.5*(ra + rp);      % GTO semi-major axis [km]

% -------------------- Event timing derived from decision vars -------------
T1   = eta1 * TOF;             % Earth → GA duration
tGA  = t0 + T1;                % GA epoch
tArr = t0 + TOF;               % Final arrival epoch
T2   = (1 - eta1) * TOF;       % GA → arrival duration
tGA2DSM   = eta2 * T2;         % GA → DSM segment duration
tDSM2Ast  = (1 - eta2) * T2;   % DSM → asteroid segment duration

% -------------------- Ephemerides: Earth at departure and GA --------------
% J2000 states wrt Sun at t0 and tGA (position [km], velocity [km/s])
D_earth1 = cspice_spkezr('Earth', t0,  'J2000','NONE','Sun');
D_earth2 = cspice_spkezr('Earth', tGA, 'J2000','NONE','Sun');
rE1 = D_earth1(1:3); vE1 = D_earth1(4:6);
rE2 = D_earth2(1:3); vE2 = D_earth2(4:6);

% -------------------- Lambert arc: Earth(t0) → Earth(tGA) -----------------
% Heliocentric Lambert to get the required departure/arrival heliocentric v.
% This delivers v1 (at rE1) and v2 (at rE2); multiple-rev branches may exist.
[v1, v2] = lambert_solver(rE1, rE2, T1, muS);

% If multiple revolution branches exist, pick the one with the smallest
% heliocentric |v| as a pragmatic tie-breaker (keeps the lowest-energy arc).
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


% -------------------- Departure Δv from an Earth GTO ----------------------
% Hyperbolic excess wrt Earth at t0:
v_inf1    = v1(:) - vE1(:);                 
v_inf1_n  = norm(v_inf1);

% Injection cost: tangential burn at GTO perigee to achieve required v∞.
% Δv_inj = v_peri,hyper - v_peri,ellipse
% (launch-vehicle Δv and parking-orbit mechanics represented succinctly).
delta_v_inj = sqrt(v_inf1_n^2 + 2*muE/rp) - sqrt(muE*(2/rp - 1/a));  

% -------------------- Unpowered GA at Earth -------------------------------
% Inbound hyperbolic excess at GA epoch:
v_inf2 = v2(:) - vE2(:);                      % inbound v∞ (km/s)
v_out  = GA_Earth_calculations(Rperi, theta_b, v_inf2, vE2, tGA);   % heliocentric velocity after GA

% -------------------- Propagate GA → DSM (Sun two-body) -------------------
% Propagation from GA epoch to the DSM epoch in heliocentric 2BP.
stDSM = cspice_prop2b(muS, [rE2; v_out(:)], tGA2DSM);
r_pre_DSM = stDSM(1:3);
v_pre_DSM = stDSM(4:6);

% -------------------- DSM → Asteroid Lambert arc --------------------------
% Target asteroid (here: 54509621 = 2024 YR4) at arrival:
D_ast = cspice_spkezr('54509621', tArr, 'J2000','NONE','Sun');
rAst  = D_ast(1:3); vAst = D_ast(4:6);

% Solve Lambert from DSM point to asteroid at tArr for the post-DSM v.
[v_post_DSM, v_sc_ast] = lambert_solver(r_pre_DSM, rAst, tDSM2Ast, muS);

% If multiple solutions exist, pick the lowest-|v| branch (as above).
if size(v_post_DSM, 1) > 1
    v0 = norm(v_post_DSM(1, :));
    v11 = norm(v_post_DSM(2, :));
    v22 = norm(v_post_DSM(3, :));
    if v0 < v11 && v0 < v22
    v_post_DSM = v_post_DSM(1, :);
    v_sc_ast = v_sc_ast(1, :);

    elseif v11 < v0 && v11 < v22
    v_post_DSM = v_post_DSM(2, :);
    v_sc_ast = v_sc_ast(2, :);

    else 
    v_post_DSM = v_post_DSM(3, :);
    v_sc_ast = v_sc_ast(3, :);
    end
end

% DSM impulse magnitude:
delta_v_DSM   = v_post_DSM(:) - v_pre_DSM(:);    % vector Δv at DSM (km/s)
delta_v_DSM_n = norm(delta_v_DSM);

% -------------------- Arrival relative speed (for flyby) -------
% Hyperbolic excess wrt the asteroid at arrival (km/s):
v_inf_arr = v_sc_ast(:) - vAst(:);
v_inf_arr_n = norm(v_inf_arr);

% -------------------- Scalar constraint thresholds ------------------------
% Minimal durations (mission phasing / numerical robustness):
minTOF   = 1*86400;        % overall minimum TOF [s]
Tmin1    = 1*86400;         % min Earth→GA duration [s]
Tmin2a   = 1*86400;         % min GA→DSM duration [s]
Tmin2b   = 1*86400;         % min DSM→arrival duration [s]

% Cap on departure v∞ to reflect launcher/escape energy constraints:
vinfMax0 = 3.873;            % [km/s]

% --- Arrival geometry  ------------------------------------------
dir_rAst = -rAst/norm(rAst); % anti-solar radial direction at arrival
dir_v_arr = -v_inf_arr/v_inf_arr_n; % opposite to arrival v∞ (approach direction)

% Calculate cosine of angle between approach direction and anti-solar direction
cos_arr_angle = dot(dir_v_arr, dir_rAst);

% Clamp to valid range for numerical stability
cos_arr_angle = max(-1, min(1, cos_arr_angle));
% -------------------- Assemble inequality constraints c <= 0 --------------
c = [];

% (1) Arrival approach geometry: prefer sunward approach (≤60° from −r̂)
% Constraint: cos(arrival_angle) ≥ cos(60°)
% Rearranged: cos(60°) - cos(arrival_angle) ≤ 0
c(end + 1) = cosd(60) - cos_arr_angle;

% (2) Overall minimum TOF.
c(end + 1) = minTOF - TOF;                 % overall TOF ≥ minTOF

% (3) Segment minimums for robustness and GA/DSM feasibility.
c(end + 1) = Tmin1  - T1;                  % Earth→GA ≥ Tmin1
c(end + 1) = Tmin2a - eta2*T2;             % GA→DSM ≥ Tmin2a
c(end + 1) = Tmin2b - (1-eta2)*T2;         % DSM→arrival ≥ Tmin2b

% (4) Departure energy cap: ||v∞,dep|| ≤ vinfMax0.
c(end + 1) = v_inf1_n - vinfMax0;          % cap departure v_inf

% (5) (Optional) Arrival v∞ cap (here 7 km/s).
c(end + 1) = v_inf_arr_n - 7;

% -------------------- No equality constraints -----------------------------
ceq = [];



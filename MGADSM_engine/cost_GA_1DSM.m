function [J] = cost_GA_1DSM(x)
% =========================================================================
% Cost function: single Earth gravity assist (unpowered) + one DSM
%
% Decision vector:
%   x = [ t0, TOF, eta1, eta2, theta_b, Rperi ]
%     t0      : departure epoch (ET, TDB seconds)
%     TOF     : total time of flight (s)
%     eta1    : fraction of TOF at which the GA occurs (Earth→GA leg = eta1*TOF)
%     eta2    : fraction along post-GA leg where DSM occurs (GA→DSM = eta2*(1-eta1)*TOF)
%     theta_b : GA B-plane turning/aiming parameter (rad) used by GA_calculations
%     Rperi   : Earth flyby periapsis radius (km)
%
% Modelling assumptions:
%   • Patched-conics, heliocentric two-body propagation between events.
%   • Planet states from SPICE (Sun/J2000, NO aberration).
%   • GA unpowered 
%   • DSM is impulsive at the chosen epoch.
% 
%
% Objective:
%   J = Δv_injection_from_GTO + Δv_DSM              (no Δv at GA; arrival Δv not charged here)

t0 = x(1); TOF = x(2); eta1 = x(3); eta2 = x(4); theta_b = x(5); Rperi = x(6);

% -------------------- Constants & GTO orbit info ----------
muS = 1.32712440018e11;   % Sun GM [km^3/s^2]
muE = 398600;             % Earth GM [km^3/s^2]
Re  = 6378.137;           % Earth radius [km]
rp  = Re + 250;           % GTO perigee radius [km]
ra  = Re + 22500;         % GTO apogee radius [km]
a   = 0.5*(ra + rp);      % GTO semi-major axis [km]

% -------------------- Event times from decision variables ------------------
T1   = eta1 * TOF;                  % Earth→GA leg duration
tGA  = t0 + T1;                     % GA epoch
tArr = t0 + TOF;                    % Arrival epoch
T2   = (1 - eta1) * TOF;            % GA→arrival total
tGA2DSM   = eta2 * T2;              % GA→DSM duration
tDSM2Ast  = (1 - eta2) * T2;        % DSM→arrival duration

% -------------------- Planetary states (Sun/J2000) ------------------------
D_earth1 = cspice_spkezr('Earth', t0,  'J2000','NONE','Sun');
D_earth2 = cspice_spkezr('Earth', tGA, 'J2000','NONE','Sun');
rE1 = D_earth1(1:3); vE1 = D_earth1(4:6);
rE2 = D_earth2(1:3); vE2 = D_earth2(4:6);

% -------------------- Lambert arc: Earth(t0) → Earth(tGA) -----------------
% Heliocentric velocities at endpoints (may return multiple branches)
[v1, v2] = lambert_solver(rE1, rE2, T1, muS);

% Branch selection: choose the arc with smallest |v1| (energy proxy)
if size(v1,1) > 1
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

% -------------------- Departure injection Δv from GTO ---------------------
% Required Earth-relative v∞ at departure:
v_inf1    = v1(:) - vE1(:);                    % [km/s]
v_inf1_n  = norm(v_inf1);

% Tangential perigee burn to achieve that v∞ from GTO:
% Δv0 = v_peri(hyper) − v_peri(ellipse) 
delta_v_inj = sqrt(v_inf1_n^2 + 2*muE/rp) - sqrt(muE*(2/rp - 1/a));  

% -------------------- Unpowered GA at Earth -------------------------------
% Inbound v∞ at GA epoch (heliocentric → planet-relative):
v_inf2 = v2(:) - vE2(:);                       % inbound v∞ at GA 
v_out  = GA_Earth_calculations(Rperi, theta_b, v_inf2, vE2, tGA);     % GA calculations
 
% -------------------- Propagate GA → DSM (Sun two-body) -----------------------
stDSM = cspice_prop2b(muS, [rE2; v_out(:)], tGA2DSM);
r_pre_DSM = stDSM(1:3);
v_pre_DSM = stDSM(4:6);

% -------------------- DSM → asteroid Lambert arc --------------------------
D_ast = cspice_spkezr('54509621', tArr, 'J2000','NONE','Sun');    % 2024 YR4
rAst  = D_ast(1:3); vAst = D_ast(4:6);

[v_post_DSM, v_sc_ast] = lambert_solver(r_pre_DSM, rAst, tDSM2Ast, muS);

% Calculate the velocity for each revolution and select the one with least
% energy
if size(v_post_DSM,1) > 1
    v0 = norm(v_post_DSM(1,:));
    v11 = norm(v_post_DSM(2,:));
    v22 = norm(v_post_DSM(3,:));
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

% -------------------- DSM impulse and arrival v∞ -------------
delta_v_DSM   = v_post_DSM(:) - v_pre_DSM(:);    % DSM vector Δv [km/s]
delta_v_DSM_n = norm(delta_v_DSM);

% Arrival hyperbolic excess (not charged in J, but useful for reporting)
v_inf_arr = v_sc_ast(:) - vAst(:);
v_inf_arr_n = norm(v_inf_arr);

% -------------------- Objective function ----------------------------------
% J = Δv0 (injection from GTO) + Δv_DSM. No Δv at GA; arrival Δv excluded.
J = delta_v_inj + delta_v_DSM_n;

end



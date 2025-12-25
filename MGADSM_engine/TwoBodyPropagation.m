function [v_pre_DSM, r_pre_DSM] = TwoBodyPropagation(t0, tf, r0, v0)

% =========================================================================
% TWO-BODY PROPAGATION (Sun-centered)
%
% Purpose:
%   Propagate a heliocentric state (r0, v0) at epoch t0 forward to epoch tf
%   under two-body dynamics, by:
%     1) extracting classical orbital elements (COEs) from (r0, v0),
%     2) advancing the mean anomaly (elliptic) or hyperbolic mean anomaly,
%     3) solving Kepler's equation for E (elliptic) or H (hyperbolic),
%     4) mapping back to true anomaly f and reconstructing the state at tf
%        via COE2RV.
%
% Inputs:
%   t0, tf : initial/final epochs (seconds, consistent with SPICE ET if used)
%   r0, v0 : initial Sun-centered inertial state in km and km/s (3x1 each)
%
% Outputs:
%   v_pre_DSM : propagated velocity at tf (km/s)
%   r_pre_DSM : propagated position at tf (km)
%
% Notes:
%   • Works for elliptic (e<1) and hyperbolic (e>1) cases; near-parabolic
%     (e≈1) uses a "nearly-elliptic" fallback to avoid singular formulas.
%   • Requires user-defined solvers Kepler(e,M,tol), keplerH(M,e) and
%     COE2RV([a e i RAAN ω θ], μ).
%   • Central body is the Sun (μ set below); units are km, s, km/s.
%   • Angles are in radians unless otherwise stated.
% =========================================================================

mu = 1.32712440018e11;                % km^3/s^2 (Sun)
TOF = tf - t0;

% helpers
clamp = @(x) max(-1, min(1, x));      % numerically safe acos/asin argument
wrap2pi = @(x) mod(x, 2*pi);          % wrap angle to [0, 2π)

% --- Orbital elements from state ------------------------------------------
r = norm(r0(:)); v = norm(v0(:));
H = cross(r0,v0); h = norm(H);     % specific angular momentum
k = [0 0 1];                       % inertial z-axis

% Inclination i
i    = acos( clamp(H(3)/h) );          % radians

n_v  = cross(k, H);  n = norm(n_v);    % Node vector and RAAN Ω
RAAN = acos( clamp(n_v(1)/n) );
if n_v(2) < 0
    RAAN = 2*pi - RAAN; 
end

% Eccentricity vector e_v and magnitude e
e_v = (1/mu) * ((v^2 - mu/r)*r0 - (dot(r0,v0)*v0));
e   = norm(e_v);

% Argument of periapsis ω
argP = acos( clamp(dot(n_v,e_v)/(n*e)) );
if e_v(3) < 0
    argP = 2*pi - argP;
end

% Semi-major axis a (works for elliptic a>0 and hyperbolic a<0)
a = 1/((2/r) - (v^2/mu));              % km

% ---- True anomaly at t0, f0 ----------------------------------------------
theta0 = acos( clamp(dot(e_v, r0 )/(e*norm(r0 ))) );
if dot(r0, v0) < 0
    theta0 = 2*pi - theta0;
end

% --- Advance anomaly to tf (solve Kepler in appropriate regime) -----------
eps_par = 1e-12;           % small buffer around e=1
f0 = theta0;               % true anomaly at t0

if e < 1 - eps_par
    % -------- Elliptic case --------
    % Eccentric anomaly at t0
    cosE0 = (e + cos(f0)) / (1 + e*cos(f0));
    sinE0 = sqrt(1 - e^2) * sin(f0) / (1 + e*cos(f0));
    E0    = atan2(sinE0, cosE0);
    E0    = wrap2pi(E0);

    % Mean motion and mean anomaly
    nmean = sqrt(mu/a^3);
    M0    = E0 - e*sinE0;
    M     = M0 + nmean * TOF;

    % Solve Kepler for E (use your existing Kepler solver)
    E = Kepler(e, M, 1e-12);

    % Back to true anomaly at t = t0+TOF
    cosf = (cos(E) - e) / (1 - e*cos(E));
    sinf = (sqrt(1 - e^2) * sin(E)) / (1 - e*cos(E));

elseif e > 1 + eps_par
    % -------- Hyperbolic case --------
    % Hyperbolic anomaly at t0
    sinhH0 = (sqrt(e^2 - 1) * sin(f0)) / (1 + e*cos(f0));
    H0     = asinh(sinhH0);

    % Hyperbolic "mean motion" and mean anomaly
    N  = sqrt(mu/(-a)^3);
    M0 = e*sinhH0 - H0;
    M  = M0 + N * TOF;

    % Solve Kepler (hyperbolic) for H
    H = keplerH(M, e);

    % Back to true anomaly
    cosf = (cosh(H) - e) / (1 - e*cosh(H));
    sinf = -(sqrt(e^2 - 1) * sinh(H)) / (1 - e*cosh(H));

else
    % -------- Near-parabolic fallback --------
    % Treat as "nearly-elliptic" to avoid singular formulas
    e_eff = max(0, 1 - 1e-12);
    cosE0 = (e_eff + cos(f0)) / (1 + e_eff*cos(f0));
    sinE0 = sqrt(1 - e_eff^2) * sin(f0) / (1 + e_eff*cos(f0));
    E0    = atan2(sinE0, cosE0);
    nmean = sqrt(mu/a^3);
    M0    = E0 - e_eff*sinE0;
    M     = M0 + nmean * TOF;
    E     = Kepler(e_eff, M, 1e-12);
    cosf  = (cos(E) - e_eff) / (1 - e_eff*cos(E));
    sinf  = (sqrt(1 - e_eff^2) * sin(E)) / (1 - e_eff*cos(E));
end

theta = atan2(sinf, cosf);       % true anomaly at tf

% --- Reconstruct state at tf from COEs ------------------------------------
coe = [a; e; i; RAAN; argP; theta];     % [a e i Ω ω θ]

[X] = COE2RV(coe, mu);                  % user routine: COE → [r; v]

r_pre_DSM = X(1:3, :);
v_pre_DSM = X(4:6, :);

end
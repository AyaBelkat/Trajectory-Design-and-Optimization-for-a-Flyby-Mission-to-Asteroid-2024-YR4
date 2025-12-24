function [v_out] = GA_calculations(Rperi, theta_b, v_inf_in, vE, et)

% GA_CALCULATIONS   Computes the outgoing heliocentric velocity after a gravity assist around Earth.
%
% INPUTS:
%   Rperi    [km]   Periapsis radius of the flyby relative to planet center
%   theta_b  [rad]  B-plane angle
%   v_inf_in [3x1]  Incoming hyperbolic excess velocity in J2000 frame
%   vE       [3x1]  Planet heliocentric velocity in J2000 frame
%   et       [s]    Ephemeris time (TDB seconds past J2000) for frame transform
%
% OUTPUT:
%   v_out    [3x1]  Outgoing heliocentric velocity in J2000 frame

muE = 398600;    % Earth's gravitational parameter [km^3/s^2]

% Magnitude of incoming excess velocity
v_inf = norm(v_inf_in);

% Transform planet pole from body-fixed to J2000 at epoch DCM

M_eq = cspice_pxform(['IAU_' 'EARTH'], 'J2000', et);   % body-fixed, equatorial -> J2000

% Define B-plane reference axes
S = v_inf_in/norm(v_inf_in);    % Incoming asymptote direction
K_eq = M_eq * [0; 0; 1];        % Planet pole in J2000
K_eq = K_eq/norm(K_eq);
T = cross(S, K_eq)/norm(cross(S, K_eq));
R = cross(S, T)/norm(cross(S, T));

% Turning angle from periapsis radius and v_inf
B = sqrt(Rperi^2 + (2*muE*Rperi)/v_inf^2);

delta = 2 * atan2(muE, B * v_inf^2);    % Turning angle

% Direction of B vector from b_plane angle
B_dir = cos(theta_b) * T + sin(theta_b) * R;
B_hat = B_dir;

% Rotate incoming v_inf by delta around B_hat (Rodrigues' formula)
v_inf_out = v_inf_in * cos(delta) + cross(B_hat, v_inf_in) * sin(delta) + B_hat * (dot(B_hat, v_inf_in)) * (1 - cos(delta));

% Convert back to heliocentric frame
v_out = vE + v_inf_out;
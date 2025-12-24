function plot_lambert_trajectory(R1, R2, V, N_rev)
% Plot the Lambert arc (elliptic or hyperbolic) between R1 and R2
% with Earth and asteroid orbits using SPICE
% Inputs:
%   R1 [km] : initial position
%   R2 [km] : final position
%   V  [km/s]: velocity at R1 (used for quadrant / flight-direction)
%
% Requires:
%   - Kepler(e, M, tol) & keplerH(e, M, tol): solves Kepler's eq (returns E for e<1, H for e>1)
%   - COE2RV(coe, mu): converts [a;e;i;RAAN;argP;theta] (rad) to r,v
%   - SPICE/MICE toolkit


mu = 1.32712440018e11;                % km^3/s^2 (Sun)
AU = 1.495978707e8;                   % 1 AU in km

% helpers
clamp   = @(x) max(-1, min(1, x));
wrap2pi = @(x) mod(x, 2*pi);

r = norm(R1); v = norm(V);
Hvec = cross(R1,V); h = norm(Hvec);
k = [0 0 1];

% Inclination, RAAN
i    = acos( clamp(Hvec(3)/h) );          % [rad]
n_v  = cross(k, Hvec);  n = norm(n_v);
RAAN = acos( clamp(n_v(1)/n) );
if n_v(2) < 0
    RAAN = 2*pi - RAAN;
end

% Eccentricity vector & e, semi-major axis a (signed)
e_v = (1/mu) * ((v^2 - mu/r)*R1 - (dot(R1,V)*V));
e   = norm(e_v);
a   = 1/((2/r) - (v^2/mu));               % km (negative for hyperbolic)

% Argument of periapsis
argP = acos( clamp(dot(n_v,e_v)/(n*e)) );
if e_v(3) < 0
    argP = 2*pi - argP;
end

% --- angles MUST stay in radians for COE2RV ---
% RAAN = RAAN * 180/pi; i = i * 180/pi; argP = argP * 180/pi;  % leave commented

% ---- Find true anomalies at R1 and R2 (zero-rev, short-way) ----
theta1 = acos( clamp(dot(e_v, R1)/(e*norm(R1))) );
if dot(R1, V) < 0
    theta1 = 2*pi - theta1;
end

theta2 = acos( clamp(dot(e_v, R2)/(e*norm(R2))) );
% choose the quadrant of theta2 consistent with flight direction (H right-hand)
if dot(cross(R2, V), Hvec) < 0
    theta2 = 2*pi - theta2;
end

% Add the revolution count (you need to pass this from your Lambert solver)
theta2 = theta2 + N_rev * 2*pi;  % where N_rev = 0, 1, 2, ...
% -------------------- Time of flight & propagation --------------------
ELL_EPS = 1e-10;
tol = 1e-10;

if e < 1 - ELL_EPS
    % -------------------- Elliptic branch --------------------
    % Eccentric anomalies
    E1 = 2*atan2( sqrt(1-e)*sin(theta1/2), sqrt(1+e)*cos(theta1/2) ); E1 = wrap2pi(E1);
    if N_rev == 0
        E2 = 2*atan2( sqrt(1-e)*sin(theta2/2), sqrt(1+e)*cos(theta2/2) ); E2 = wrap2pi(E2);
    else
        E2 = 2*atan2( sqrt(1-e)*sin(theta2/2), sqrt(1+e)*cos(theta2/2) );
    end

    % Mean anomalies & mean motion
    M1 = E1 - e*sin(E1);
    M2 = E2 - e*sin(E2);
    nmean = sqrt(mu/a^3);

    % Short-way delta M (0..2pi)
    if N_rev == 0
        dM = wrap2pi(M2 - M1);
    else
        dM = M2 - M1;
    end
    
    TOF = dM / nmean;                      % seconds

    % Propagate from M1 to M2
    tspan = linspace(0, TOF, 1000);
    X = zeros(6, numel(tspan));
    for tt = 1:numel(tspan)
        M = M1 + nmean*tspan(tt);
        E = Kepler(e, M, tol);                         % elliptic solution
        theta = 2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );
        coe = [a; e; i; RAAN; argP; theta];
        X(:,tt) = COE2RV(coe, mu);
    end
    orbit_type = 'Elliptic';
    P = 2*pi*sqrt(a^3/mu);   % period
    %fprintf('%s orbit. Period = %.6f s, TOF = %.6f s\n', orbit_type, P, TOF);

elseif e > 1 + ELL_EPS
    % -------------------- Hyperbolic branch --------------------
    % Hyperbolic eccentric anomalies from true anomaly
    % tanh(H/2) = sqrt((e-1)/(e+1)) * tan(theta/2)
    H1 = 2*atanh( sqrt((e-1)/(e+1)) * tan(theta1/2) );
    H2 = 2*atanh( sqrt((e-1)/(e+1)) * tan(theta2/2) );

    % Hyperbolic mean anomalies & "mean motion"
    M1 = e*sinh(H1) - H1;
    M2 = e*sinh(H2) - H2;
    nh = sqrt(mu/(-a)^3);                   % note: a<0 so -a>0

    dM = M2 - M1;                           % should be >0 along flight
    if dM <= 0
        % If the chosen quadrants produced a negative TOF, flip sign
        warning('Computed negative hyperbolic TOF; taking absolute value.');
        dM = abs(dM);
    end
    TOF = dM / nh;

    % Propagate from M1 to M2
    tspan = linspace(0, TOF, 1000);
    X = zeros(6, numel(tspan));
    for tt = 1:numel(tspan)
        M = M1 + nh*tspan(tt);
        Hh = KeplerH(e, M, tol);                            % hyperbolic solution
        % true anomaly from H:
        % tan(theta/2) = sqrt((e+1)/(e-1)) * tanh(H/2)
        th = 2*atan( sqrt((e+1)/(e-1)) * tanh(Hh/2) );
        % ensure th is in (-pi,pi) (no wrap for hyperbola)
        if th > pi, th = th - 2*pi;
        elseif th < -pi, th = th + 2*pi;
        end
        coe = [a; e; i; RAAN; argP; th];
        X(:,tt) = COE2RV(coe, mu);
    end
    orbit_type = 'Hyperbolic';
    fprintf('%s orbit. TOF = %.6f s (period undefined)\n', orbit_type, TOF);

else
    % -------------------- Parabolic (near-e==1) --------------------
    error(['Parabolic case (e≈1) not handled here. ', ...
           'Use Barker''s equation or a universal-variable propagator.']);
end

% --- Clip trajectory at actual arrival point (closest to R2) --------------
d = sqrt(sum((X(1:3,:) - R2).^2, 1));   % distance to target along the arc
[~, idx] = min(d);                       % index of closest approach
X = X(:, 1:idx);                         % keep only up to that point

% ---------- Get Earth and asteroid orbits using SPICE ----------
fprintf('Calculating Earth and asteroid orbits using SPICE...\n');

% Time setup
start_date = datetime(2028, 1, 1);
t0 = cspice_str2et(datestr(start_date));
duration_days_earth = 365.25;                % One Earth year
duration_days_asteroid = 1457.59;    % Lambert TOF + buffer

tf1 = t0 + duration_days_earth*86400;
tf2 = t0 + duration_days_asteroid*86400;

% Create time vectors
N_samples = 1000;
time_vector_E = linspace(t0, tf1, N_samples);
time_vector_A = linspace(t0, tf2, N_samples);

% Get Earth orbit
X_earth = zeros(6, N_samples);
for i = 1:N_samples
    earth_state = cspice_spkezr('Earth', time_vector_E(i), 'J2000', 'NONE', 'Sun');
    X_earth(:, i) = earth_state;
end

% Get 2024 YR4 orbit
X_yr4 = zeros(6, N_samples);
for i = 1:N_samples
    try
        yr4_state = cspice_spkezr('54509621', time_vector_A(i), 'J2000', 'NONE', 'Sun');
        X_yr4(:, i) = yr4_state;
    catch
        X_yr4(:, i) = [0;0;0;0;0;0];  % Zero if unavailable
    end
end

% ---- Plot ----
R_s = 696340;   % Sun radius [km]

[Xe,Ye,Ze] = sphere(5000);
figure('units','normalized','outerposition',[0 0 1 1])
surf(R_s*Xe, R_s*Ye, R_s*Ze, 'EdgeColor','none','FaceColor','#F80'); hold on; axis equal
xlabel('X, SCI (km)'); ylabel('Y, SCI (km)'); zlabel('Z, SCI (km)');
grid on

% Plot Earth orbit
plot3(X_earth(1,:), X_earth(2,:), X_earth(3,:), 'b-', 'LineWidth', 1.5);

% Plot asteroid orbit (if valid data exists)
valid_idx = ~all(X_yr4(1:3,:) == 0, 1);
if sum(valid_idx) > 1
    plot3(X_yr4(1,valid_idx), X_yr4(2,valid_idx), X_yr4(3,valid_idx), 'r-', 'LineWidth', 1.5);
end

% Plot Lambert arc
plot3(X(1,:), X(2,:), X(3,:), 'k', 'LineWidth', 2);
P1 = plot3(X(1,1),  X(2,1),  X(3,1),  'ok', 'MarkerFaceColor','b', 'DisplayName','Departure: Earth');   % start at R1
P2 = plot3(R2(1),   R2(2),   R2(3),   'ok', 'MarkerFaceColor','r', 'DisplayName', 'Arrival: Asteroid 2024 YR4');   % end at R2

% Add Earth and asteroid orbit labels
P3 = plot3(X_earth(1,1), X_earth(2,1), X_earth(3,1), 'MarkerFaceColor','b', 'DisplayName', 'Earth Orbit');
if sum(valid_idx) > 1
    P4 = plot3(X_yr4(1,find(valid_idx,1)), X_yr4(2,find(valid_idx,1)), X_yr4(3,find(valid_idx,1)), 'MarkerFaceColor','r', 'DisplayName', '2024 YR4 Orbit');
    legend([P1(1),P2(1),P3(1),P4(1)],'Location','best')
else
    legend([P1(1),P2(1),P3(1)],'Location','best')
end

title(sprintf('%s Lambert Arc with Planetary Orbits', orbit_type));
hold off

% Clean up SPICE (optional)
% cspice_kclear;

end

% === Earth → 2024 YR4 Lambert scan (C3/TOF pork-chop inputs) ===============
% High-level:
%   - Builds departure/arrival date grids
%   - Pulls heliocentric J2000 states with SPICE (Earth, 2024 YR4)
%   - Solves Lambert for each (dep, arr) pair → possibly multiple solutions
%   - Records C3 = ||v∞,dep||^2, TOF [days], and arrival/departure v∞ magnitudes
%   - Filters/NaNs large values for plotting; calls user plotting utilities
%
% Frames/units/aberration:
%   - Frame: 'J2000' (ICRF), center: 'Sun', ABCORR: 'NONE'
%   - Distances [km], time [s], velocities [km/s]

clc
clear
format longG

mu = 1.32712440018e11;    % Sun gravitational parameter [km^3/s^2]

% --- SPICE setup: load MICE paths and required kernels ---------------------
addpath("C:\Users\wbook\Desktop\Dissertation files\Dissertation_codes\mice\mice\src\mice");
addpath("C:\Users\wbook\Desktop\Dissertation files\Dissertation_codes\mice\mice\lib");

cspice_furnsh('naif0012.tls.pc');               % Leap seconds (time conversions)
cspice_furnsh('54509621_2024yr4_data.bsp');     % SPK for asteroid 2024 YR4 (NAIF ID 54509621)

cspice_furnsh('de430.bsp');                     % Planetary ephemerides 

% --- Date grids (UTC calendar) ---------------------------------------------
departure_dates = datetime(2031,7,1):days(1):datetime(2032,4,15);
arrival_dates = datetime(2032,6,1):days(1):datetime(2032,12,15);

numDepartures = length(departure_dates);
numArrivals = length(arrival_dates);

% --- Convert calendar dates → SPICE ephemeris time (ET, TDB seconds) -------
for i = 1:numDepartures
    depEt(i) = cspice_str2et(datestr(departure_dates(i)));
end

for j = 1:numArrivals
    arrEt(j) = cspice_str2et(datestr(arrival_dates(j)));
end

% Prealloc Arrays will grow as used below.
C3  = NaN(numArrivals, numDepartures);
TOF = NaN(numArrivals, numDepartures);

% --- Main Lambert scan ------------------------------------------------------
for i = 1:numDepartures
    for j = 1:numArrivals

        % Skip non-causal pairs (arrival must be after departure)
        if arrEt(j) - depEt(i) <= 0
            continue;
        end

        % Heliocentric states at requested epochs (target frame: Sun/J2000)
        % D_earth, D_ast: [6x1] (pos; vel) in km and km/s at ET
        D_earth = cspice_spkezr('Earth', depEt(i), 'J2000', 'NONE', 'Sun');
        D_ast = cspice_spkezr('54509621',  arrEt(j), 'J2000', 'NONE', 'Sun');
        tof = arrEt(j) - depEt(i);     % Time of flight [s]

        % Extract position/velocity endpoints for Lambert
        r1 = D_earth(1:3);
        r2 = D_ast(1:3);

        vE = D_earth(4:6);
        vA = D_ast(4:6);
        
        % Solve Lambert (may return multiple branches: 0-rev, multi-rev)
        [v1, v2] = lambert_solver(r1, r2, tof, mu);    % Lambert solver
        
        % --- Record "first" solution (index 1) --------------------------------
        % Departure hyperbolic excess: v∞,dep = v1 − v_Earth(dep)
        v_inf_1(j, i, 1) = norm(v1(1,:) - vE');
        C3(j, i, 1) = v_inf_1(j, i, 1)^2;      % C3 = |v∞,dep|^2  [km^2/s^2]

        TOF(j, i, 1) = tof/86400;              % Time of flight [days]
        vv1(j, i, 1) = norm(v1(1,:));          % |v1| (diagnostic)
        vv2(j, i, 1) = norm(v2(1,:));          % |v2| (diagnostic)

        % Arrival hyperbolic excess: v∞,arr = v_Ast(arr) − v2
        v_inf_2(j, i, 1) = norm(vA' - v2(1,:));

        % --- If multiple Lambert solutions exist, store up to 3 ----------------
        if size(v1, 1) > 1
            % Norms of additional branches (diagnostics)
            vv1(j, i, 2) = norm(v1(2,:));
            vv2(j, i, 2) = norm(v2(2,:));
            vv1(j, i, 3) = norm(v1(3,:));
            vv2(j, i, 3) = norm(v2(3,:));

            % Corresponding v∞ at departure/arrival for branches 2 and 3
            v_inf_1(j, i, 2) = norm(v1(2,:) - vE');
            v_inf_1(j, i, 3) = norm(v1(3,:) - vE');

            v_inf_2(j, i, 2) = norm(v2(2,:) - vA');
            v_inf_2(j, i, 3) = norm(v2(3,:) - vA');

            % C3 and TOF mirrors for these branches (same TOF grid)
            C3(j, i, 2) = v_inf_1(j, i, 2)^2;
            C3(j, i, 3) = v_inf_1(j, i, 3)^2;

            TOF(j, i, 2) = tof/86400;
            TOF(j, i, 3) = tof/86400;
        end
    end
end

% --- Plot-prep: mask large C3 and keep TOF consistent where C3 is NaN -----
C3_plot = C3;
C3_plot(C3_plot >= 11) = NaN;    % Visual cap: C3 ≥ 11 km^2/s^2 hidden

tof_plot = TOF;                  % Copy TOF and NaN-mask wherever C3 is NaN
for i = 1:size(C3_plot, 1)
    for j = 1:size(C3_plot, 2)
        for k = 1:size(C3_plot, 3)
            if isnan(C3_plot(i,j,k))
                tof_plot(i, j, k) = NaN;
            end
        end
    end
end

% Contour levels (tune for readability)
levels = 0:2:10;          % C3 contours [km^2/s^2]
levels_tof = 0:20:700;    % TOF contours [days]

% --- User plotting utility: C3 vs. (dep, arr) with TOF overlays ------------
plotC3(C3_plot, tof_plot, levels, levels_tof, departure_dates, arrival_dates)


% --- Delta-v summaries from v∞ magnitudes ----------------
% Note: delta_v here := |v∞,dep| + |v∞,arr|
delta_v = v_inf_1 + v_inf_2;

% Visual caps for plots (hide extreme tails)
delta_v_plot = delta_v;
delta_v_plot(delta_v_plot >= 21) = NaN;

delta_v1_plot = v_inf_1;
delta_v1_plot(delta_v1_plot >= 11) = NaN;

delta_v2_plot = v_inf_2;
delta_v2_plot(delta_v2_plot >= 11) = NaN;

Total_delta_v(departure_dates, arrival_dates, delta_v_plot, C3)

delta_v_2_plot(delta_v2_plot, departure_dates, arrival_dates, C3)

delta_v_1_plot(delta_v1_plot, departure_dates, arrival_dates, C3)


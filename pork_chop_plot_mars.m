clc
clear



mu = 1.32712440018e11;   % Sun gravitational constant [km^3/s^2]

% === SPICE/MICE setup =====================================================
% Adds MICE toolboxes to the MATLAB path and loads core kernels:
%  - naif0012.tls.pc : leap seconds (UTC↔TDB conversions)
%  - de430.bsp       : JPL planetary/lunar ephemerides
addpath("C:\Users\wbook\Desktop\Dissertation files\Dissertation_codes\mice\mice\src\mice");
addpath("C:\Users\wbook\Desktop\Dissertation files\Dissertation_codes\mice\mice\lib");

cspice_furnsh('naif0012.tls.pc');    % Time system support
cspice_furnsh('de430.bsp');          % Ephemerides for Earth/Mars (center: Sun)

% === Date grids (UTC calendar) ============================================
% Earth→Mars verification case (2005–2007):
%   - departures every 2 days
%   - arrivals every 5 days
departure_dates = datetime(2005,6,20,0,0,0):days(2):datetime(2005,11,7,0,0,0);
arrival_dates = datetime(2005,12,1,0,0,0):days(5):datetime(2007,2,24,0,0,0);

numDepartures = length(departure_dates);
numArrivals = length(arrival_dates);

% Convert calendar datetimes → SPICE ephemeris time (ET, TDB seconds)
for i = 1:numDepartures
    depEt(i) = cspice_str2et(datestr(departure_dates(i)));
end

for j = 1:numArrivals
    arrEt(j) = cspice_str2et(datestr(arrival_dates(j)));
end

% === Allocate pork-chop arrays ============================================
% C3  : launch energy (km^2/s^2) based on departure v-infinity
% TOF : time of flight (days)
C3  = NaN(numArrivals, numDepartures, 3);
TOF = NaN(numArrivals, numDepartures, 3);

v_inf_1 = NaN(numArrivals, numDepartures, 3);
v_inf_2 = NaN(numArrivals, numDepartures, 3);
vv1     = NaN(numArrivals, numDepartures, 3);
vv2     = NaN(numArrivals, numDepartures, 3);

% === Main Lambert scan loop ==============================================
% For each (departure, arrival) pair:
%   1) skip non-causal pairs (TOF ≤ 0)
%   2) pull heliocentric states (Sun/J2000, no aberration)
%   3) solve Lambert
%   4) compute departure/arrival v∞ and derived metrics
for i = 1:numDepartures
    for j = 1:numArrivals

        % Skip non-causal pairs (arrival must be after departure)
        if arrEt(j) - depEt(i) <= 0
            continue;
        end

        % States wrt Sun in J2000 (ICRF), km and km/s, at given ET
        D_earth = cspice_spkezr('Earth', depEt(i), 'J2000', 'NONE', 'Sun');
        D_mars = cspice_spkezr('4',  arrEt(j), 'J2000', 'NONE', 'Sun');     % NAIF ID 4 = Mars

        tof = arrEt(j) - depEt(i);    % TOF [s]
        
        % Heliocentric states at requested epochs (target frame: Sun/J2000)
        % D_earth, D_ast: [6x1] (pos; vel) in km and km/s at ET
        r1 = D_earth(1:3);
        r2 = D_mars(1:3);
        
        vE = D_earth(4:6);
        vM = D_mars(4:6);
        
        % Solve Lambert (may return multiple branches: 0-rev, multi-rev)
        [v1, v2] = lambert_solver(r1, r2, tof, mu);     % Lambert solver
        
        % --- Record "first" solution (index 1) --------------------------------
        % Departure hyperbolic excess: v∞,dep = v1 − v_Earth(dep)
        v_inf_1(j, i, 1) = norm(v1(1,:) - vE');
        C3(j, i, 1) = v_inf_1(j, i, 1)^2;        % Launch energy C3 = |v∞,dep|^2

        TOF(j, i, 1) = tof/86400;        % TOF [days] (for plotting/contours)
        vv1(j, i, 1) = norm(v1(1,:));    % |v1| (diagnostic)
        vv2(j, i, 1) = norm(v2(1,:));    % |v2| (diagnostic)

        % Arrival hyperbolic excess: v∞,arr = v_Mars(arr) − v2
        v_inf_2(j, i, 1) = norm(vM' - v2(1,:));
        
        % --- If multiple Lambert solutions exist, store up to 3 ------------
        if size(v1, 1) > 1

            % Norms of additional branches (diagnostics)
            vv1(j, i, 2) = norm(v1(2,:));
            vv2(j, i, 2) = norm(v2(2,:));
            vv1(j, i, 3) = norm(v1(3,:));
            vv2(j, i, 3) = norm(v2(3,:));

            % Corresponding v∞ at departure/arrival for branches 2 and 3
            v_inf_1(j, i, 2) = norm((v1(2,:)) - vE');
            v_inf_1(j, i, 3) = norm((v1(3,:)) - vE');

            v_inf_2(j, i, 2) = norm(v2(2,:) - vM');
            v_inf_2(j, i, 3) = norm(v2(3,:) - vM');

            % C3 and TOF mirrors for these branches (same TOF grid)
            C3(j, i, 2) = v_inf_1(j, i, 2)^2;
            C3(j, i, 3) = v_inf_1(j, i, 3)^2;

            TOF(j, i, 2) = tof/86400;
            TOF(j, i, 3) = tof/86400;
        end

    end
end


% === Plot-prep masks for readability ======================================
% Cap the visualised C3 domain to < 31 km^2/s^2; carry NaNs into TOF plot
C3_plot = C3;
C3_plot(C3_plot >= 31) = NaN;

TOF_plot = TOF;
for i = 1:size(C3_plot, 1)
    for j = 1:size(C3_plot, 2)
        for k = 1:size(C3_plot, 3)
            if isnan(C3_plot(i,j,k))
                TOF_plot(i, j, k) = NaN;
            end
        end
    end
end

% === Pork-chop plot (C3 with TOF overlays) =================================
figure;
% Build the meshgrid of datenums for contour axes
[X, Y] = meshgrid(datenum(departure_dates), datenum(arrival_dates));

% Contour levels (tune for readability)
levels = 0:2:30;         % C3 [km^2/s^2]
levels_tof = 0:25:600;   % TOF [days]

% Draw contours per branch (k), overlay TOF in magenta
for k = 1:size(C3, 3)
    contour(X, Y, C3_plot(:,:,k), levels, 'LineWidth',1, 'ShowText','on');
    hold on;
    contour(X, Y, TOF_plot(:,:,k), levels_tof,'m-', 'LineWidth',3, 'ShowText','on');
end
hold off

% Colouring & scale
colormap(jet);
colorbar;
clim([10 30])     % Clamp visible colour range (here focusing on typical E→M values)

% Axes limits & ticks in date format (UTC calendar on both axes)
ax = gca;

ax.XLim = datenum([departure_dates(1), departure_dates(end)]);
ax.YLim = datenum([arrival_dates(1),   arrival_dates(end)]);

ax.XTick = linspace(ax.XLim(1), ax.XLim(2), 20);
ax.YTick = linspace(ax.YLim(1), ax.YLim(2), 20);

datetick(ax, 'x', 'yyyy-mm-dd','keeplimits', 'keepticks');
datetick(ax, 'y', 'yyyy-mm-dd', 'keeplimits','keepticks');

ax.XTickLabelRotation = 45;
% ax.FontSize = 9;

set(ax, 'XGrid','on','YGrid','on');

xlabel('Departure Dates');
ylabel('Arrival Dates');
title('Earth–Mars Pork Chop: C_3 < 30 km^2/s^2');


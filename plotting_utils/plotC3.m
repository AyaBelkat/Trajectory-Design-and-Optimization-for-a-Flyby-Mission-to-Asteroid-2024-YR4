function plotC3(C3_plot, levels, departure_dates, arrival_dates)
% === Function: plotC3 =====================================================
% Purpose:
%   Generate a pork chop plot to visualise launch energy requirements (C3) 
%   for transfers between Earth and asteroid 2024 YR4.
%
% Inputs:
%   C3_plot        : 3D matrix of launch energy (km^2/s^2) solutions from Lambert solver
%   TOF            : Time of flight (days), used optionally for overlay contours
%   levels         : Contour levels for C3 (e.g. 0–10 km^2/s^2)
%   levels_tof     : Contour levels for TOF (days) [currently commented]
%   departure_dates: Array of departure epochs (datetime format)
%   arrival_dates  : Array of arrival epochs (datetime format)
%
% Outputs:
%   None (generates a 2D contour figure)
% ==========================================================================

% --- Meshgrid for departure/arrival dates (datenum converts datetime -> serial) ---
[X, Y] = meshgrid(datenum(departure_dates), datenum(arrival_dates));

% --- Plot C3 contours for each revolution branch (k = index) ---------------
figure;
for k = 1:size(C3_plot,3)
    contour(X, Y, C3_plot(:,:,k),levels,'LineWidth',1);
    hold on
end
hold off

% --- Colour map and scale -------------------------------------------------
colormap(jet);
colormap(jet);
colorbar;
clim([0 10]);    % Limit C3 colour scale between 0–10 km^2/s^2

% --- Axes formatting ------------------------------------------------------
ax = gca;
ax.XLim = datenum([departure_dates(1), departure_dates(end)]);
ax.YLim = datenum([arrival_dates(1),   arrival_dates(end)]);

% Define ~20 ticks for readability on both axes
ax.XTick = linspace(ax.XLim(1), ax.XLim(2), 20);
ax.YTick = linspace(ax.YLim(1), ax.YLim(2), 20);

% Format axes as calendar dates
datetick(ax, 'x', 'yyyy-mm-dd','keeplimits', 'keepticks');
datetick(ax, 'y', 'yyyy-mm-dd', 'keeplimits','keepticks');


ax.XTickLabelRotation = 45;   % Rotate labels for clarity
ax.FontSize = 9;

set(ax, 'XGrid','on','YGrid','on');

% --- Labels and title -----------------------------------------------------
xlabel('Departure Dates');
ylabel('Arrival Dates');
title('Earth–2024YR4 Pork Chop: C_3 < 10 km^2/s^2');

% --- Enable interactive click on plot (retrieve departure/arrival pair) ---
ax = gca;
set(allchild(ax), 'HitTest','off', 'PickableParts','none');                 % children ignore clicks
set(ax, 'HitTest','on', 'PickableParts','all', 'ButtonDownFcn', @ClickOnContour);  % axes handles clicks



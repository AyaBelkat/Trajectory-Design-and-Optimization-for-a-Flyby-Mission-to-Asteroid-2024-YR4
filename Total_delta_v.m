function Total_delta_v(departure_dates, arrival_dates, delta_v, C3)

figure;

[X, Y] = meshgrid(datenum(departure_dates), datenum(arrival_dates));

%levels = 0:1:10;

for k = 1:size(C3,3)
    contour(X, Y, delta_v(:,:,k), 'LineWidth',1, 'ShowText','on');
    hold on;
end
hold off

colormap(jet);
colorbar;
%clim([0 10])      % clamp the colorbar to [0,25]

ax = gca;

ax.XLim = datenum([departure_dates(1), departure_dates(end)]);
ax.YLim = datenum([arrival_dates(1),   arrival_dates(end)]);

ax.XTick = linspace(ax.XLim(1), ax.XLim(2), 20);
ax.YTick = linspace(ax.YLim(1), ax.YLim(2), 20);


datetick(ax, 'x', 'yyyy-mm-dd','keeplimits', 'keepticks');
datetick(ax, 'y', 'yyyy-mm-dd', 'keeplimits','keepticks');


ax.XTickLabelRotation = 45;
ax.FontSize = 9;

set(ax, 'XGrid','on','YGrid','on');

xlabel('Departure Dates');
ylabel('Arrival Dates');
title('Earth–2024YR4 total delta v plot');
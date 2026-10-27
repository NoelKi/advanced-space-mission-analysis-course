function plotOrbit(r_eci, utc, titleStr)
% PLOTORBIT Zeichnet die Bahn in 3D (ECI) und die Bodenspur nebeneinander.
%
% Eingaben:
%   r_eci     Positionsvektoren (3xN) im Inertialsystem [km]
%   utc       Zeitpunkte (1xN datetime, UTC)
%   titleStr  Titel der Abbildung
    c = earthConstants();
    figure('Name', titleStr, 'Position', [100 100 1300 520]);
    tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
    title(tl, titleStr);

    % 3D-Bahn mit Erdkugel
    nexttile;
    [X, Y, Z] = sphere(60);
    surf(c.Re*X, c.Re*Y, c.Re*Z, 'FaceColor', [0.55 0.75 1], ...
        'EdgeColor', [0.4 0.55 0.8], 'EdgeAlpha', 0.3, 'FaceAlpha', 0.6);
    hold on;
    plot3(r_eci(1,:), r_eci(2,:), r_eci(3,:), 'r', 'LineWidth', 0.8);
    plot3(r_eci(1,1), r_eci(2,1), r_eci(3,1), 'ko', 'MarkerFaceColor', 'g');
    axis equal; grid on; view(135, 25);
    xlabel('x_{ICRF} [km]'); ylabel('y_{ICRF} [km]'); zlabel('z_{ICRF} [km]');
    title('Bahn im Inertialsystem');
    legend('Erde', 'Bahn', 'Start', 'Location', 'northeast');

    % Bodenspur; an Sprüngen über den Datumsgrenzmeridian aufgetrennt
    [lat, lon] = eci2geodetic(r_eci, utc);
    gx = geoaxes(tl);
    gx.Layout.Tile = 2;
    geobasemap(gx, 'darkwater');
    hold(gx, 'on');
    breaks = [0, find(abs(diff(lon)) > 180), numel(lon)];
    for s = 1:numel(breaks)-1
        idx = breaks(s)+1 : breaks(s+1);
        geoplot(gx, lat(idx), lon(idx), 'r', 'LineWidth', 0.8);
    end
    geoplot(gx, lat(1), lon(1), 'ko', 'MarkerFaceColor', 'g');
    geolimits(gx, [-90 90], [-180 180]);
    title(gx, 'Bodenspur');
end

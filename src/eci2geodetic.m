function [lat, lon, h] = eci2geodetic(r_eci, utc)
% ECI2GEODETIC Rechnet Inertialpositionen in geodätische Koordinaten um.
%
% Die Erdrotation wird über die Greenwich Mean Sidereal Time (GMST, IAU 1982)
% berücksichtigt; Präzession, Nutation und Polbewegung werden vernachlässigt.
% Die geodätische Breite wird iterativ auf dem WGS84-Ellipsoid bestimmt.
%
% Eingaben:
%   r_eci  Positionsvektoren (3xN) im Inertialsystem [km]
%   utc    Zeitpunkte (1xN datetime, UTC)
%
% Rückgabewerte:
%   lat  Geodätische Breite [deg] (1xN)
%   lon  Geographische Länge [deg] in [-180, 180] (1xN)
%   h    Höhe über dem Ellipsoid [km] (1xN)
    c = earthConstants();
    d = juliandate(utc) - 2451545.0;                 % Tage seit J2000
    T = d / 36525;                                   % Julianische Jahrhunderte
    gmst = deg2rad(280.46061837 + 360.98564736629*d + 0.000387933*T.^2 - T.^3/38710000);

    % Drehung um die z-Achse: ECI -> ECEF
    x =  cos(gmst).*r_eci(1,:) + sin(gmst).*r_eci(2,:);
    y = -sin(gmst).*r_eci(1,:) + cos(gmst).*r_eci(2,:);
    z =  r_eci(3,:);

    lon = rad2deg(atan2(y, x));

    e2  = c.f * (2 - c.f);                           % Exzentrizität^2 des Ellipsoids
    rxy = hypot(x, y);
    phi = atan2(z, rxy * (1 - e2));                  % Startwert
    for k = 1:5
        Nphi = c.Re ./ sqrt(1 - e2*sin(phi).^2);     % Querkrümmungsradius
        h    = rxy ./ cos(phi) - Nphi;
        phi  = atan2(z, rxy .* (1 - e2*Nphi./(Nphi + h)));
    end
    lat = rad2deg(phi);
end

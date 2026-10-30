classdef Earth
    properties (Constant)
        mu      = 398600.4418               % Gravitationsparameter der Erde [km^3/s^2]
        Re      = 6378.137                  % Äquatorradius (WGS84) [km]
        f       = 1/298.257223563           % Abplattung (WGS84) [-]
        J2      = 1.08262668e-3             % Zonaler Koeffizient J2 [-]
        wE      = 7.2921159e-5              % Rotationsrate der Erde [rad/s]
        raanSSO = 2*pi / (365.2422 * 86400) % Knotendrift einer SSO, 360° pro trop. Jahr [rad/s]
    end
end
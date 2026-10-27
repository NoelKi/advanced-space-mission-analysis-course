function c = earthConstants()
% EARTHCONSTANTS Liefert die in allen Aufgaben verwendeten Erdkonstanten.
%
% Rückgabewert:
%   c.mu      Gravitationsparameter der Erde [km^3/s^2]
%   c.Re      Äquatorradius (WGS84) [km]
%   c.f       Abplattung (WGS84) [-]
%   c.J2      Zonaler Koeffizient J2 [-]
%   c.wE      Rotationsrate der Erde [rad/s]
%   c.raanSSO Erforderliche Knotendrift einer SSO, 360° pro tropischem Jahr [rad/s]
    c.mu      = 398600.4418;
    c.Re      = 6378.137;
    c.f       = 1/298.257223563;
    c.J2      = 1.08262668e-3;
    c.wE      = 7.2921159e-5;
    c.raanSSO = 2*pi / (365.2422 * 86400);
end

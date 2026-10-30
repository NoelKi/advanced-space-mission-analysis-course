classdef Moon
    properties (Constant)
        mu      = 4902.800066               % Gravitationsparameter des Mondes (DE430) [km^3/s^2]
        Re      = 1738.0                    % Referenzradius der GRAIL-Schwerefeldmodelle [km]
        f       = 0.0012                    % Abplattung (1738.1 km äquatorial, 1736.0 km polar) [-]
        J2      = 2.0322e-4                 % Zonaler Koeffizient J2 (GRAIL GRGM1200A) [-]
        wE      = 2*pi / (27.321661 * 86400) % Rotationsrate des Mondes, gebunden: 1 Umdrehung pro sid. Monat [rad/s]
        raanSSO = 2*pi / (365.2422 * 86400) % Knotendrift einer sonnensynchronen Bahn, 360° pro trop. Jahr [rad/s]
    end
end

function RAANdot = raanRateJ2(a, e, i)
% RAANRATEJ2 Säkulare Knotendrift durch die J2-Störung.
%
%   dRAAN/dt = -3/2 * n * J2 * (Re/p)^2 * cos(i)
%
% Eingaben:
%   a  Große Halbachse [km]
%   e  Exzentrizität [-]
%   i  Inklination [rad]
%
% Rückgabewert:
%   RAANdot  Drift der Rektaszension des aufsteigenden Knotens [rad/s]
    c = earthConstants();
    n = sqrt(c.mu / a^3);              % mittlere Bewegung
    p = a * (1 - e^2);                 % Semi-latus rectum
    RAANdot = -1.5 * n * c.J2 * (c.Re/p)^2 * cos(i);
end

function i = ssoInclination(a, e)
% SSOINCLINATION Inklination einer sonnensynchronen Bahn.
%
% Die J2-Knotendrift (siehe raanRateJ2) wird gleich der mittleren
% Bewegung der Sonne gesetzt (360° pro tropischem Jahr) und nach cos(i)
% aufgelöst:
%   cos(i) = -RAANdot_SSO / (3/2 * n * J2 * (Re/p)^2)
%
% Eingaben:
%   a  Große Halbachse [km]
%   e  Exzentrizität [-]
%
% Rückgabewert:
%   i  Inklination [rad] (retrograd, > 90°)
    c = earthConstants();
    cosi = c.raanSSO / raanRateJ2(a, e, 0);   % raanRateJ2(a,e,0) = -3/2*n*J2*(Re/p)^2
    if abs(cosi) > 1
        error('ssoInclination:noSolution', ...
            'Für a = %.1f km und e = %.3f existiert keine SSO.', a, e);
    end
    i = acos(cosi);
end

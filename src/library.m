% Library by Kieran Noel Mai for 
% (advanced) space mission analysis 

% Beispiel: a [km], e [-], Winkel in rad
kep = [7000, 0.01, deg2rad(98), deg2rad(30), deg2rad(45), deg2rad(10)];
[r_eci, v_eci] = kep2cart(kep)

% Beispiel: Kepler-Gleichung M = E - e*sin(E) nach E lösen
e = 0.01;
M = deg2rad(10);
[E, iter] = newton(@(E) E - e*sin(E) - M, @(E) 1 - e*cos(E), M)

function [r_eci, v_eci] = kep2cart(kep, mu)
% KEP2CART Rechnet Kepler-Elemente in kartesische Zustandsvektoren um.
%
% Rückgabewerte:
%   r_eci  Positionsvektor (3x1) im Inertialsystem (ECI), gemessen vom
%          Mittelpunkt des Zentralkörpers [km]
%   v_eci  Geschwindigkeitsvektor (3x1) im Inertialsystem (ECI) [km/s]
%
% Die Einheiten gelten für a in km und mu in km^3/s^2.
    arguments
        kep (1,6) double {mustBeReal, mustBeFinite}
        mu  (1,1) double {mustBePositive} = 398600.4418   % km^3/s^2 (Erde)
    end
    a     = kep(1);
    e     = kep(2);
    i     = kep(3);
    RAAN  = kep(4);
    omega = kep(5);
    nu    = kep(6);
    
    p = a * (1 - e^2);                 % Semi-latus rectum
    r = p / (1 + e*cos(nu));           % Bahnradius
    
    r_pf = r * [cos(nu); sin(nu); 0];
    v_pf = sqrt(mu/p) * [-sin(nu); e + cos(nu); 0];
    
    R3_RAAN  = [ cos(RAAN)  sin(RAAN) 0; -sin(RAAN)  cos(RAAN) 0; 0 0 1];
    R1_i     = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
    R3_omega = [ cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];
    
    Q = (R3_omega * R1_i * R3_RAAN)';  % perifokal -> inertial
    
    r_eci = Q * r_pf;
    v_eci = Q * v_pf;
end

function [x, iter] = newton(f, df, x0, opts)
% NEWTON Sucht eine Nullstelle von f mit dem Newton-Verfahren.
%
% Eingaben:
%   f             Funktions-Handle der Funktion, deren Nullstelle gesucht wird
%   df            Funktions-Handle der Ableitung f'
%   x0            Startwert
%   opts.Tol      Abbruch, sobald |x_(n+1) - x_n| < Tol (Standard: 1e-12)
%   opts.MaxIter  Maximale Anzahl an Iterationen (Standard: 50)
%
% Rückgabewerte:
%   x     Gefundene Nullstelle
%   iter  Anzahl der benötigten Iterationen
    arguments
        f  (1,1) function_handle
        df (1,1) function_handle
        x0 (1,1) double {mustBeReal, mustBeFinite}
        opts.Tol     (1,1) double {mustBePositive} = 1e-12
        opts.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 50
    end
    x = x0;
    for iter = 1:opts.MaxIter
        dfx = df(x);
        if dfx == 0
            error('newton:zeroDerivative', 'Ableitung ist null bei x = %g.', x);
        end
        dx = f(x) / dfx;                   % Newton-Schritt: x_(n+1) = x_n - f/f'
        x  = x - dx;
        if abs(dx) < opts.Tol
            return
        end
    end
    error('newton:noConvergence', 'Keine Konvergenz nach %d Iterationen.', opts.MaxIter);
end
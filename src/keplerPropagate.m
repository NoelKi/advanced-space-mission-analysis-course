function [r_eci, v_eci] = keplerPropagate(kep0, t, opts)
% KEPLERPROPAGATE Propagiert eine Keplerbahn analytisch über die Zeit.
%
% Für jeden Zeitpunkt wird die mittlere Anomalie fortgeschrieben, die
% Kepler-Gleichung M = E - e*sin(E) mit dem Newton-Verfahren nach E gelöst
% und daraus die wahre Anomalie und der Zustandsvektor berechnet.
% Optional wird die säkulare Knotendrift durch J2 berücksichtigt.
%
% Eingaben:
%   kep0     Kepler-Elemente zur Epoche t0: [a e i RAAN omega nu], a in km,
%            Winkel in rad
%   t        Zeitpunkte seit t0 [s] (1xN)
%   opts.J2  true: RAAN driftet mit raanRateJ2 (Standard: false)
%
% Rückgabewerte:
%   r_eci  Positionsvektoren (3xN) im Inertialsystem [km]
%   v_eci  Geschwindigkeitsvektoren (3xN) im Inertialsystem [km/s]
    arguments
        kep0 (1,6) double {mustBeReal, mustBeFinite}
        t    (1,:) double {mustBeReal, mustBeFinite}
        opts.J2 (1,1) logical = false
    end
    c = earthConstants();
    a     = kep0(1);
    e     = kep0(2);
    i     = kep0(3);
    RAAN0 = kep0(4);
    omega = kep0(5);
    nu0   = kep0(6);

    n  = sqrt(c.mu / a^3);                                     % mittlere Bewegung
    E0 = 2*atan2(sqrt(1-e)*sin(nu0/2), sqrt(1+e)*cos(nu0/2));  % wahre -> exz. Anomalie
    M0 = E0 - e*sin(E0);

    RAANdot = 0;
    if opts.J2
        RAANdot = raanRateJ2(a, e, i);
    end

    N = numel(t);
    r_eci = zeros(3, N);
    v_eci = zeros(3, N);
    for k = 1:N
        M  = mod(M0 + n*t(k), 2*pi);
        E  = newton(@(E) E - e*sin(E) - M, @(E) 1 - e*cos(E), M);
        nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2)); % exz. -> wahre Anomalie
        RAAN = RAAN0 + RAANdot*t(k);
        [r_eci(:,k), v_eci(:,k)] = kep2cart([a e i RAAN omega nu], c.mu);
    end
end

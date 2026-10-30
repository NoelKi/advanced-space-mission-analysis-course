%% Space Mission Analysis WS 2026/27 – Assignment 1
% Kieran Noel Mai
%
% Annahmen laut Aufgabenstellung: omega = 0°, RAAN0 = 0°, nu0 = 0°.
% Startepoche t0 frei gewählt. Alle Funktionen liegen in der Klasse OrbitKit.
clear; close all; clc;

t0 = datetime(2026, 12, 24, 0, 0, 0, 'TimeZone', 'UTC');   % Startepoche (UTC)

%% Aufgabe 1a – Kepler-Propagation
% Implementiert in OrbitKit.keplerPropagate (Kepler-Gleichung mit newton gelöst).
% Zeitschritt: 10 s ergeben ca. 580 Punkte pro Umlauf. Da die Lösung analytisch
% ist, beeinflusst der Schritt nur die Auflösung der Plots, nicht die Genauigkeit.
dt = 10;   % [s]

%% Aufgabe 1b – SSO-Inklination und ungestörte Propagation
hp = 550;                          % Perigäumshöhe [km]
e  = 0.01;
a  = (Earth.Re + hp) / (1 - e);  % aus r_p = a*(1 - e)
i  = OrbitKit.ssoInclination(a, e);
T  = 2*pi * sqrt(a^3 / Earth.mu); % Umlaufzeit [s]

kep0 = [a, e, i, 0, 0, 0];         % [a e i RAAN omega nu]
nRev = 30;
t    = 0:dt:nRev*T;
utc  = t0 + seconds(t);

fprintf('Große Halbachse a      = %.3f km\n', a);
fprintf('SSO-Inklination i      = %.4f°\n', rad2deg(i));
fprintf('Umlaufzeit T           = %.2f min\n', T/60);
fprintf('Propagationsdauer      = %.2f h (%d Punkte)\n', t(end)/3600, numel(t));

[r_kep, v_kep] = OrbitKit.keplerPropagate(kep0, t);   % Positionen (ICRF) in r_kep [km]

%% Aufgabe 1c – Plots ungestörte Bahn
OrbitKit.plotOrbit(r_kep, utc, 'Aufgabe 1c: Ungestörte Keplerbahn (SSO, 30 Umläufe)');

%% Aufgabe 1d – Kepler-Propagation mit J2-Knotendrift
% OrbitKit.keplerPropagate(..., J2=true) lässt RAAN mit raanRateJ2 linear driften.
RAANdot = OrbitKit.raanRateJ2(a, e, i);
fprintf('\nKnotendrift durch J2    = %.4f °/Tag\n', rad2deg(RAANdot)*86400);
fprintf('Soll-Drift für SSO      = %.4f °/Tag\n', rad2deg(Earth.raanSSO)*86400);

%% Aufgabe 1e – Propagation der SSO mit J2
[r_kepJ2, v_kepJ2] = OrbitKit.keplerPropagate(kep0, t, J2=true);   % Positionen (ICRF) in r_kepJ2 [km]

dr = vecnorm(r_kepJ2 - r_kep);
fprintf('RAAN nach %d Umläufen   = %.4f°\n', nRev, rad2deg(RAANdot*t(end)));
fprintf('Max. Positionsabweichung J2 vs. ungestört = %.1f km\n', max(dr));

%% Aufgabe 1f – Plots Bahn mit J2-Knotendrift
OrbitKit.plotOrbit(r_kepJ2, utc, 'Aufgabe 1f: Keplerbahn mit J2-Knotendrift (SSO, 30 Umläufe)');

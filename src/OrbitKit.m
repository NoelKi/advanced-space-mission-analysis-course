classdef OrbitKit
% ORBITKIT Sammlung aller Funktionen für (Advanced) Space Mission Analysis.
% Kieran Noel Mai
%
% Alle Funktionen sind statische Methoden und werden ohne Objekt aufgerufen,
% z. B. OrbitKit.kep2cart(kep). Die Erdkonstanten liegen in der Klasse
% Earth, z. B. Earth.mu.


    methods (Static)

        %% Numerik

        function [x, iter] = newton(f, df, x0, opts)
        % NEWTON Sucht eine Nullstelle von f mit dem Newton-Verfahren.
        %
        % Eingaben:
        %   f             Funktions-Handle der Funktion, deren Nullstelle gesucht wird
        %   df            Funktions-Handle der Ableitung f'
        %   x0            Startwert
        %   opts.Tol      Abbruch, sobald |x_(n+1) - x_n| < Tol (Tolerenz) (Standard: 1e-12)
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

        %% Koordinatentransformationen

        function [r_eci, v_eci] = kep2cart(kep, mu)
        % KEP2CART Rechnet Kepler-Elemente in kartesische Zustandsvektoren um.
        %
        % Eingaben:
        %   kep  [a e i RAAN omega nu], a in km, Winkel in rad
        %   mu   Gravitationsparameter [km^3/s^2] (Standard: Erde)
        %
        % Rückgabewerte:
        %   r_eci  Positionsvektor (3x1) im Inertialsystem (ECI), gemessen vom
        %          Mittelpunkt des Zentralkörpers [km]
        %   v_eci  Geschwindigkeitsvektor (3x1) im Inertialsystem (ECI) [km/s]
            arguments
                kep (1,6) double {mustBeReal, mustBeFinite}
                mu  (1,1) double {mustBePositive} = Earth.mu
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
            d = juliandate(utc) - 2451545.0;                 % Tage seit J2000
            T = d / 36525;                                   % Julianische Jahrhunderte
            gmst = deg2rad(280.46061837 + 360.98564736629*d + 0.000387933*T.^2 - T.^3/38710000);

            % Drehung um die z-Achse: ECI -> ECEF
            x =  cos(gmst).*r_eci(1,:) + sin(gmst).*r_eci(2,:);
            y = -sin(gmst).*r_eci(1,:) + cos(gmst).*r_eci(2,:);
            z =  r_eci(3,:);

            lon = rad2deg(atan2(y, x));

            Re  = Earth.Re;
            e2  = Earth.f * (2 - Earth.f);   % Exzentrizität^2 des Ellipsoids
            rxy = hypot(x, y);
            phi = atan2(z, rxy * (1 - e2));                  % Startwert
            for k = 1:5
                Nphi = Re ./ sqrt(1 - e2*sin(phi).^2);       % Querkrümmungsradius
                h    = rxy ./ cos(phi) - Nphi;
                phi  = atan2(z, rxy .* (1 - e2*Nphi./(Nphi + h)));
            end
            lat = rad2deg(phi);
        end

        %% Bahnmechanik

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
            n = sqrt(Earth.mu / a^3);   % mittlere Bewegung
            p = a * (1 - e^2);                     % Semi-latus rectum
            RAANdot = -1.5 * n * Earth.J2 * (Earth.Re/p)^2 * cos(i);
        end

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
            % raanRateJ2(a,e,0) = -3/2*n*J2*(Re/p)^2
            cosi = Earth.raanSSO / OrbitKit.raanRateJ2(a, e, 0);
            if abs(cosi) > 1
                error('ssoInclination:noSolution', ...
                    'For a = %.1f km and e = %.3f there is no SSO solution.', a, e);
            end
            i = acos(cosi);
        end

        function dr_rtn = eci2rtn(dr_eci, r_ref, v_ref)
        % ECI2RTN Zerlegt Differenzvektoren in das RTN-System einer Referenzbahn.
        %
        % R radial (Richtung r), N normal zur Bahnebene (Richtung r x v),
        % T entlang der Bahn (N x R).
        %
        % Eingaben:
        %   dr_eci  Differenzvektoren (3xN) im Inertialsystem
        %   r_ref   Positionen der Referenzbahn (3xN) [km]
        %   v_ref   Geschwindigkeiten der Referenzbahn (3xN) [km/s]
        %
        % Rückgabewert:
        %   dr_rtn  Komponenten [R; T; N] (3xN), gleiche Einheit wie dr_eci
            eR = r_ref ./ vecnorm(r_ref);
            eN = cross(r_ref, v_ref);
            eN = eN ./ vecnorm(eN);
            eT = cross(eN, eR);
            dr_rtn = [dot(dr_eci, eR); dot(dr_eci, eT); dot(dr_eci, eN)];
        end

        %% Bewegungsgleichungen

        function dy = eomTwoBody(~, y)
        % EOMTWOBODY Ungestörte Bewegungsgleichung (Zweikörperproblem).
        %
        %   r'' = -mu/|r|^3 * r
        %
        % Eingaben:
        %   y  Zustand [r; v] (6x1) im Inertialsystem [km, km/s]
        %
        % Rückgabewert:
        %   dy  Zeitableitung [v; a] (6x1) [km/s, km/s^2]
            r  = y(1:3);
            dy = [y(4:6); -Earth.mu / norm(r)^3 * r];
        end

        function dy = eomJ2(~, y)
        % EOMJ2 Bewegungsgleichung mit J2-Störung (Abplattung der Erde).
        %
        %   r'' = -mu/|r|^3 * r + a_J2,  mit
        %   a_J2 = -3/2 * J2 * mu * Re^2 / |r|^5 * [x*(1 - 5 z^2/r^2)
        %                                           y*(1 - 5 z^2/r^2)
        %                                           z*(3 - 5 z^2/r^2)]
        % Die z-Achse des Inertialsystems wird als Rotationsachse der Erde
        % angenommen (Präzession und Nutation vernachlässigt).
        %
        % Eingaben:
        %   y  Zustand [r; v] (6x1) im Inertialsystem [km, km/s]
        %
        % Rückgabewert:
        %   dy  Zeitableitung [v; a] (6x1) [km/s, km/s^2]
            mu  = Earth.mu;
            r   = y(1:3);
            rn  = norm(r);
            z2  = (r(3) / rn)^2;
            kJ2 = -1.5 * Earth.J2 * mu * Earth.Re^2 / rn^5;
            aJ2 = kJ2 * [r(1)*(1 - 5*z2); r(2)*(1 - 5*z2); r(3)*(3 - 5*z2)];
            dy  = [y(4:6); -mu / rn^3 * r + aJ2];
        end

        %% Propagation

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
            mu    = Earth.mu;
            a     = kep0(1);
            e     = kep0(2);
            i     = kep0(3);
            RAAN0 = kep0(4);
            omega = kep0(5);
            nu0   = kep0(6);

            n  = sqrt(mu / a^3);                                       % mittlere Bewegung
            E0 = 2*atan2(sqrt(1-e)*sin(nu0/2), sqrt(1+e)*cos(nu0/2));  % wahre -> exz. Anomalie
            M0 = E0 - e*sin(E0);

            RAANdot = 0;
            if opts.J2
                RAANdot = OrbitKit.raanRateJ2(a, e, i);
            end

            N = numel(t);
            r_eci = zeros(3, N);
            v_eci = zeros(3, N);
            for k = 1:N
                M  = mod(M0 + n*t(k), 2*pi);
                E  = OrbitKit.newton(@(E) E - e*sin(E) - M, @(E) 1 - e*cos(E), M);
                nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2)); % exz. -> wahre Anomalie
                RAAN = RAAN0 + RAANdot*t(k);
                [r_eci(:,k), v_eci(:,k)] = OrbitKit.kep2cart([a e i RAAN omega nu], mu);
            end
        end

        %% Darstellung

        function plotOrbit(r_eci, utc, titleStr)
        % PLOTORBIT Zeichnet die Bahn in 3D (ECI) und die Bodenspur nebeneinander.
        %
        % Eingaben:
        %   r_eci     Positionsvektoren (3xN) im Inertialsystem [km]
        %   utc       Zeitpunkte (1xN datetime, UTC)
        %   titleStr  Titel der Abbildung
            Re = Earth.Re;
            figure('Name', titleStr, 'Position', [100 100 1300 520]);
            tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
            title(tl, titleStr);

            % 3D-Bahn mit Erdkugel
            nexttile;
            [X, Y, Z] = sphere(60);
            surf(Re*X, Re*Y, Re*Z, 'FaceColor', [0.55 0.75 1], ...
                'EdgeColor', [0.4 0.55 0.8], 'EdgeAlpha', 0.3, 'FaceAlpha', 0.6);
            hold on;
            plot3(r_eci(1,:), r_eci(2,:), r_eci(3,:), 'color', [1 0 1], 'LineWidth', 0.8);
            plot3(r_eci(1,1), r_eci(2,1), r_eci(3,1), 'ko', 'MarkerFaceColor', 'r');
            axis equal; grid on; view(135, 25);
            xlabel('x_{ICRF} [km]'); ylabel('y_{ICRF} [km]'); zlabel('z_{ICRF} [km]');
            title('Bahn im Inertialsystem');
            legend('Earth', 'Orbit', 'Start', 'Location', 'northeast');

            % Bodenspur; an Sprüngen über den Datumsgrenzmeridian aufgetrennt
            [lat, lon] = OrbitKit.eci2geodetic(r_eci, utc);
            gx = geoaxes(tl);
            gx.Layout.Tile = 2;
            geobasemap(gx, 'satellite');
            hold(gx, 'on');
            breaks = [0, find(abs(diff(lon)) > 180), numel(lon)];
            for s = 1:numel(breaks)-1
                idx = breaks(s)+1 : breaks(s+1);
                geoplot(gx, lat(idx), lon(idx), 'color', [1 0 1], 'LineWidth', 0.8);
            end
            geoplot(gx, lat(1), lon(1), 'ko', 'MarkerFaceColor', 'r');
            geolimits(gx, [-90 90], [-180 180]);
            title(gx, 'Bodenspur');
        end

    end
end

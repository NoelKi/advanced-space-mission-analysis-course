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
nRev = 30;                         % Orbit anzahl
t    = 0:dt:nRev*T;                % 
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

%% Aufgabe 2a – Runge-Kutta-Löser
% Implementiert in der Klasse RungeKutta (Butcher-Tableau als Eingabe,
% Schrittweitensteuerung für eingebettete Verfahren). Die Bewegungsgleichungen
% stehen in OrbitKit.eomTwoBody und OrbitKit.eomJ2.
rk45  = RungeKutta.rkf45();
rk78  = RungeKutta.rkf78();
rkOpt = {'RelTol', 1e-10, 'AbsTol', 1e-10, 'Step', 1};    % gemeinsame Einstellungen (Startschritt 1 s laut Vorlesung)

[r0, v0] = OrbitKit.kep2cart(kep0);
y0 = [r0; v0];

%% Aufgabe 2b – Ungestörte Bewegungsgleichung (RKF7(8))
[t_rk, Y_rk] = rk78.solve(@OrbitKit.eomTwoBody, [0 nRev*T], y0, rkOpt{:});
r_rk = Y_rk(1:3,:);                                      % Positionen (ICRF) [km]

% Vergleich mit der analytischen Keplerbahn aus 1b an denselben Zeitpunkten
dr_2b = vecnorm(r_rk - OrbitKit.keplerPropagate(kep0, t_rk));
fprintf('\n2b: %d Schritte, max. Abweichung zur Keplerbahn = %.2f m\n', ...
    numel(t_rk) - 1, 1e3*max(dr_2b));

%% Aufgabe 2c – Bewegungsgleichung mit J2 (RKF7(8))
[t_rkJ2, Y_rkJ2] = rk78.solve(@OrbitKit.eomJ2, [0 nRev*T], y0, rkOpt{:});
r_rkJ2 = Y_rkJ2(1:3,:);                                  % Positionen (ICRF) [km]

% Vergleich mit der Kepler-Propagation mit Knotendrift aus 1e
dr_2c = vecnorm(r_rkJ2 - OrbitKit.keplerPropagate(kep0, t_rkJ2, J2=true));
fprintf('2c: %d Schritte, max. Abweichung zu Kepler mit J2-Knotendrift = %.1f km\n', ...
    numel(t_rkJ2) - 1, max(dr_2c));

figure('Name', 'Aufgabe 2b/2c');
semilogy(t_rk/3600, 1e3*dr_2b, t_rkJ2/3600, 1e3*dr_2c);
grid on; xlabel('t [h]'); ylabel('|\Deltar| [m]');
legend('2b: RKF7(8) ungestört – Kepler (1b)', '2c: RKF7(8) mit J2 – Kepler mit J2-Knotendrift (1e)', ...
    'Location', 'southeast');
title('Aufgabe 2: Abweichung der numerischen Lösung von der Kepler-Propagation');

%% Aufgabe 3 – Vergleich der Propagationsverfahren
% Annahmen: "Höhe" der SSO = Perigäumshöhe (wie in Aufgabe 1). Für die GTO ist
% keine Inklination vorgegeben; gewählt wird i = 6° (typische Ariane-GTO ab Kourou),
% Apogäum auf GEO-Radius. Alle Verfahren mit J2, Dauer 36 h.
rGEO = (Earth.mu / Earth.wE^2)^(1/3);                    % Radius der geostationären Bahn [km]
rp3  = Earth.Re + 650;                                   % Perigäumsradius beider Bahnen [km]
aSSO = rp3 / (1 - 0.01);
aGTO = (rp3 + rGEO) / 2;
eGTO = (rGEO - rp3) / (rGEO + rp3);

orbits = struct( ...
    'name', {'SSO', 'GTO'}, ...
    'kep',  {[aSSO, 0.01, OrbitKit.ssoInclination(aSSO, 0.01), 0, 0, 0], ...
             [aGTO, eGTO, deg2rad(6), 0, 0, 0]});

tEnd3   = 36 * 3600;                 % Propagationsdauer [s]
tq      = 0:60:tEnd3;                % gemeinsames Vergleichsgitter (Kepler und RK) [s]
nRuns   = 5;                         % Wiederholungen für die Laufzeitmessung (Median)
methodNames = {'Kepler (J2)', 'RKF4(5)', 'RKF7(8)'};
solvers = {rk45, rk78};
nOrb    = numel(orbits);

figStep = figure('Name', 'Aufgabe 3ii: Schrittweiten', 'Position', [50 50 1300 700]);
tlStep  = tiledlayout(figStep, 2, nOrb, 'TileSpacing', 'compact');
figDiff = figure('Name', 'Aufgabe 3iii: Zustandsdifferenzen', 'Position', [80 80 1300 900]);
tlDiff  = tiledlayout(figDiff, 3, nOrb, 'TileSpacing', 'compact');

for o = 1:nOrb
    kep = orbits(o).kep;
    [r0, v0] = OrbitKit.kep2cart(kep);
    y0 = [r0; v0];
    fprintf('\n3 – %s: a = %.1f km, e = %.4f, i = %.2f°, T = %.2f h\n', orbits(o).name, ...
        kep(1), kep(2), rad2deg(kep(3)), 2*pi*sqrt(kep(1)^3/Earth.mu)/3600);

    % Perturbierte Kepler-Propagation auf dem Vergleichsgitter
    tRun = zeros(1, nRuns);
    for k = 1:nRuns
        tic;
        [rK, vK] = OrbitKit.keplerPropagate(kep, tq, J2=true);
        tRun(k) = toc;
    end
    orbits(o).time = median(tRun);
    orbits(o).X    = {[rK; vK]};
    orbits(o).nSteps = numel(tq) - 1;  orbits(o).nRej = 0;  orbits(o).nFev = 0;

    % RKF4(5) und RKF7(8) mit J2-Bewegungsgleichung
    hStep = gobjects(1, 2*numel(solvers));
    for s = 1:numel(solvers)
        for k = 1:nRuns
            tic;
            [ts, Ys, st] = solvers{s}.solve(@OrbitKit.eomJ2, [0 tEnd3], y0, rkOpt{:});
            tRun(k) = toc;
        end
        orbits(o).time(end+1)   = median(tRun);
        orbits(o).X{end+1}      = solvers{s}.evaluateAt(@OrbitKit.eomJ2, ts, Ys, tq);
        orbits(o).nSteps(end+1) = st.nAccepted;
        orbits(o).nRej(end+1)   = st.nRejected;
        orbits(o).nFev(end+1)   = st.nFevals;

        % 3ii: akzeptierte und abgelehnte Schrittweiten
        ax  = nexttile(tlStep, o);
        col = ax.ColorOrder(s, :);
        hStep(2*s-1) = semilogy(ax, ts(1:end-1)/3600, st.h, '.-', 'MarkerSize', 4, 'Color', col);
        hold(ax, 'on');
        hStep(2*s) = semilogy(ax, [NaN, st.tRejected/3600], [NaN, st.hRejected], 'x', ...
            'MarkerSize', 8, 'LineWidth', 1.2, 'Color', 0.6*col);
    end
    grid(ax, 'on'); ylim(ax, 'padded'); xlabel(ax, 't [h]'); ylabel(ax, 'h [s]');
    title(ax, sprintf('%s: Schrittweiten', orbits(o).name));
    legend(hStep, 'RKF4(5) akzeptiert', 'RKF4(5) abgelehnt', 'RKF7(8) akzeptiert', ...
        'RKF7(8) abgelehnt', 'Location', 'best');

    % Bahnhöhe zur Erklärung der Schrittweiten
    ax = nexttile(tlStep, o + nOrb);
    plot(ax, tq/3600, vecnorm(rK) - Earth.Re);
    grid(ax, 'on'); xlabel(ax, 't [h]'); ylabel(ax, 'Höhe [km]');
    title(ax, sprintf('%s: Bahnhöhe', orbits(o).name));

    % 3iii: Differenzen in Position und Geschwindigkeit zwischen den Verfahren
    pairs  = {[1 3], [2 3], [1 2]};
    labels = cellfun(@(p) sprintf('%s – %s', methodNames{p(1)}, methodNames{p(2)}), pairs, ...
        'UniformOutput', false);
    X = orbits(o).X;
    for row = 1:2
        ax = nexttile(tlDiff, o + (row-1)*nOrb);
        rows = 3*(row-1) + (1:3);                        % 1: Position, 2: Geschwindigkeit
        lineStyles = {'-', '-', '--'};
        for p = 1:numel(pairs)
            semilogy(ax, tq/3600, 1e3*vecnorm(X{pairs{p}(1)}(rows,:) - X{pairs{p}(2)}(rows,:)), ...
                lineStyles{p});
            hold(ax, 'on');
        end
        grid(ax, 'on'); xlabel(ax, 't [h]');
        if row == 1
            ylabel(ax, '|\Deltar| [m]');
            title(ax, sprintf('%s: Positionsdifferenz', orbits(o).name));
            legend(ax, labels, 'Location', 'southeast');
        else
            ylabel(ax, '|\Deltav| [m/s]');
            title(ax, sprintf('%s: Geschwindigkeitsdifferenz', orbits(o).name));
        end
    end

    % Kepler – RKF7(8) im RTN-System der RKF7(8)-Bahn
    ax = nexttile(tlDiff, o + 2*nOrb);
    plot(ax, tq/3600, OrbitKit.eci2rtn(X{1}(1:3,:) - X{3}(1:3,:), X{3}(1:3,:), X{3}(4:6,:)));
    grid(ax, 'on'); xlabel(ax, 't [h]'); ylabel(ax, '\Deltar [km]');
    title(ax, sprintf('%s: Kepler (J2) – RKF7(8) im RTN-System', orbits(o).name));
    legend(ax, 'radial', 'entlang der Bahn', 'normal', 'Location', 'best');
end

%% Aufgabe 3iv – Direkter Vergleich aller drei Verfahren
fprintf('\n%-4s %-12s %12s %9s %10s %10s %18s\n', 'Bahn', 'Verfahren', ...
    'Laufzeit[ms]', 'Schritte', 'abgelehnt', 'f-Aufrufe', '|dr| zu RKF7(8)[m]');
for o = 1:nOrb
    for m = 1:numel(methodNames)
        drEnd = 1e3 * norm(orbits(o).X{m}(1:3,end) - orbits(o).X{3}(1:3,end));
        fprintf('%-4s %-12s %12.1f %9d %10d %10d %18.3g\n', orbits(o).name, methodNames{m}, ...
            1e3*orbits(o).time(m), orbits(o).nSteps(m), orbits(o).nRej(m), orbits(o).nFev(m), drEnd);
    end
end

figure('Name', 'Aufgabe 3iv: Direkter Vergleich', 'Position', [120 120 1100 420]);
tiledlayout(1, 2, 'TileSpacing', 'compact');
orbitCat = categorical({orbits.name}, {orbits.name});   % Reihenfolge wie in orbits
nexttile;
bar(orbitCat, 1e3*vertcat(orbits.time));
grid on; ylabel('Laufzeit [ms]');
legend(methodNames, 'Location', 'northwest'); title('Laufzeit (Median aus 5 Läufen)');
nexttile;
nFev = vertcat(orbits.nFev);
bar(orbitCat, nFev(:, 2:3));
grid on; ylabel('Auswertungen der Bewegungsgleichung');
legend(methodNames(2:3), 'Location', 'northwest'); title('Rechenaufwand der RK-Verfahren');

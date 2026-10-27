classdef RungeKutta
% RUNGEKUTTA Expliziter Runge-Kutta-Löser für beliebige Butcher-Tableaus.
% Kieran Noel Mai
%
% Ein Objekt beschreibt ein RK-Verfahren über sein Butcher-Tableau (A, b, c).
% Eingebettete Verfahren haben zusätzlich Gewichte bErr niedrigerer Ordnung:
% propagiert wird mit der höheren Ordnung (b), die Differenz zur niedrigeren
% Ordnung (bErr) dient als Fehlerschätzer für die Schrittweitensteuerung.
%
% Beispiel:
%   rk = RungeKutta.rkf45();
%   [t, Y, stats] = rk.solve(@(t, y) -y, [0 10], 1, RelTol=1e-8);
%
% Vordefinierte Verfahren: RungeKutta.euler, .rk4, .rkf45, .rkf78

    properties (SetAccess = immutable)
        Name     (1,:) char             % Bezeichnung des Verfahrens
        A        (:,:) double           % Koeffizientenmatrix (s x s, strikt untere Dreiecksmatrix)
        b        (1,:) double           % Gewichte der Lösung (Ordnung Order)
        c        (1,:) double           % Knoten (Zeitpunkte der Stufen)
        Order    (1,1) double           % Ordnung von b
        bErr     (1,:) double           % Gewichte für die Fehlerschätzung (leer: nicht eingebettet)
        OrderErr (1,:) double           % Ordnung von bErr
    end

    properties (Dependent)
        Stages      % Anzahl der Stufen s
        IsEmbedded  % true, wenn ein Fehlerschätzer vorhanden ist
    end

    methods

        function obj = RungeKutta(A, b, c, opts)
        % RUNGEKUTTA Erzeugt ein Verfahren aus seinem Butcher-Tableau.
        %
        % Eingaben:
        %   A              Koeffizientenmatrix (s x s)
        %   b              Gewichte (1 x s)
        %   c              Knoten (1 x s)
        %   opts.Order     Ordnung von b
        %   opts.bErr      Gewichte niedrigerer Ordnung (1 x s), optional
        %   opts.OrderErr  Ordnung von bErr, optional
        %   opts.Name      Bezeichnung, optional
            arguments
                A (:,:) double
                b (1,:) double
                c (1,:) double
                opts.Order    (1,1) double {mustBeInteger, mustBePositive}
                opts.bErr     (1,:) double = []
                opts.OrderErr (1,:) double = []
                opts.Name     (1,:) char = 'Runge-Kutta'
            end
            s = numel(b);
            if ~isequal(size(A), [s s]) || numel(c) ~= s
                error('RungeKutta:size', 'A muss %dx%d und c 1x%d sein.', s, s, s);
            end
            if any(triu(A) ~= 0, 'all')
                error('RungeKutta:implicit', 'Nur explizite Verfahren (A strikt unten) werden unterstützt.');
            end
            if ~isempty(opts.bErr) && (numel(opts.bErr) ~= s || isempty(opts.OrderErr))
                error('RungeKutta:embedded', 'bErr braucht %d Einträge und eine Angabe OrderErr.', s);
            end
            obj.A        = A;
            obj.b        = b;
            obj.c        = c;
            obj.Order    = opts.Order;
            obj.bErr     = opts.bErr;
            obj.OrderErr = opts.OrderErr;
            obj.Name     = opts.Name;
        end

        function s = get.Stages(obj)
            s = numel(obj.b);
        end

        function tf = get.IsEmbedded(obj)
            tf = ~isempty(obj.bErr);
        end

        function [t, Y, stats] = solve(obj, f, tspan, y0, opts)
        % SOLVE Integriert y' = f(t, y) von tspan(1) bis tspan(2).
        %
        % Eingebettete Verfahren arbeiten mit Schrittweitensteuerung, alle
        % anderen (oder Adaptive=false) mit fester Schrittweite opts.Step.
        %
        % Eingaben:
        %   f              Funktions-Handle f(t, y), liefert eine Spalte wie y
        %   tspan          [t0 tEnd] mit tEnd > t0
        %   y0             Anfangszustand (Spaltenvektor)
        %   opts.Step      Feste bzw. erste Schrittweite (Standard: (tEnd-t0)/1000)
        %   opts.Adaptive  Schrittweitensteuerung (Standard: IsEmbedded)
        %   opts.RelTol    Relative Toleranz (Standard: 1e-10)
        %   opts.AbsTol    Absolute Toleranz (Standard: 1e-10)
        %   opts.HMin      Kleinste zulässige Schrittweite (Standard: 1e-10)
        %   opts.HMax      Größte zulässige Schrittweite (Standard: Inf)
        %   opts.Safety    Sicherheitsfaktor der Steuerung (Standard: 0.9)
        %   opts.FacMin    Kleinster Änderungsfaktor von h (Standard: 0.2)
        %   opts.FacMax    Größter Änderungsfaktor von h (Standard: 5)
        %
        % Rückgabewerte:
        %   t      Zeitpunkte der akzeptierten Schritte (1 x N)
        %   Y      Zustände (n x N), Spalte k gehört zu t(k)
        %   stats  .h          akzeptierte Schrittweiten (1 x N-1)
        %          .tRejected  Startzeitpunkte abgelehnter Schritte
        %          .hRejected  abgelehnte Schrittweiten
        %          .nAccepted, .nRejected, .nFevals
            arguments
                obj
                f     (1,1) function_handle
                tspan (1,2) double {mustBeReal, mustBeFinite}
                y0    (:,1) double {mustBeReal, mustBeFinite}
                opts.Step     (1,1) double {mustBePositive} = (tspan(2) - tspan(1)) / 1000
                opts.Adaptive (1,1) logical = obj.IsEmbedded
                opts.RelTol   (1,1) double {mustBePositive} = 1e-10
                opts.AbsTol   (1,1) double {mustBePositive} = 1e-10
                opts.HMin     (1,1) double {mustBePositive} = 1e-10
                opts.HMax     (1,1) double {mustBePositive} = Inf
                opts.Safety   (1,1) double {mustBePositive} = 0.9
                opts.FacMin   (1,1) double {mustBePositive} = 0.2
                opts.FacMax   (1,1) double {mustBePositive} = 5
            end
            if tspan(2) <= tspan(1)
                error('RungeKutta:tspan', 'tspan(2) muss größer als tspan(1) sein.');
            end
            if opts.Adaptive && ~obj.IsEmbedded
                error('RungeKutta:notEmbedded', '%s hat keinen Fehlerschätzer.', obj.Name);
            end
            t0   = tspan(1);
            tEnd = tspan(2);

            % Ausgabe vorbelegen; wird bei Bedarf verdoppelt
            cap = 1024;
            t = zeros(1, cap);
            Y = zeros(numel(y0), cap);
            t(1) = t0;
            Y(:,1) = y0;
            N = 1;

            stats = struct('h', [], 'tRejected', [], 'hRejected', [], ...
                'nAccepted', 0, 'nRejected', 0, 'nFevals', 0);
            tRej = zeros(1, 0);
            hRej = zeros(1, 0);

            tk = t0;
            yk = y0;
            h  = min(opts.Step, opts.HMax);
            while tk < tEnd
                hk = min(h, tEnd - tk);          % letzten Schritt auf tEnd kürzen
                [yNew, err] = obj.step(f, tk, yk, hk, opts);
                stats.nFevals = stats.nFevals + obj.Stages;

                if opts.Adaptive
                    % Neue Schrittweite aus dem Fehler der eingebetteten Lösung:
                    % h_neu = h * Safety * (1/err)^(1/(q+1)), q = Ordnung des Schätzers
                    if err == 0
                        fac = opts.FacMax;
                    else
                        fac = opts.Safety * err^(-1/(obj.OrderErr + 1));
                    end
                    h = min(opts.HMax, hk * min(opts.FacMax, max(opts.FacMin, fac)));

                    if err > 1                    % Schritt ablehnen und wiederholen
                        tRej(end+1) = tk;         %#ok<AGROW>
                        hRej(end+1) = hk;         %#ok<AGROW>
                        if h < opts.HMin
                            error('RungeKutta:hMin', ...
                                'Schrittweite unter HMin (%g) bei t = %g.', opts.HMin, tk);
                        end
                        continue
                    end
                end

                % Schritt akzeptieren
                tk = tk + hk;
                yk = yNew;
                N  = N + 1;
                if N > cap
                    cap = 2*cap;
                    t(cap) = 0;
                    Y(:, cap) = 0;
                end
                t(N)   = tk;
                Y(:,N) = yk;
            end

            t = t(1:N);
            Y = Y(:, 1:N);
            stats.h         = diff(t);
            stats.tRejected = tRej;
            stats.hRejected = hRej;
            stats.nAccepted = N - 1;
            stats.nRejected = numel(hRej);
        end

    end

    methods (Access = private)

        function [yNew, err] = step(obj, f, t, y, h, opts)
        % STEP Ein RK-Schritt der Länge h; err ist der skalierte Fehler (<= 1: ok).
            s = obj.Stages;
            K = zeros(numel(y), s);
            for j = 1:s
                K(:,j) = f(t + obj.c(j)*h, y + h * K(:,1:j-1) * obj.A(j,1:j-1)');
            end
            yNew = y + h * K * obj.b';

            err = 0;
            if opts.Adaptive
                e     = h * K * (obj.b - obj.bErr)';                       % Fehlerschätzung
                scale = opts.AbsTol + opts.RelTol * max(abs(y), abs(yNew));
                err   = max(abs(e) ./ scale);
            end
        end

    end

    methods (Static)

        function obj = euler()
        % EULER Explizites Euler-Verfahren (Ordnung 1).
            obj = RungeKutta(0, 1, 0, Order=1, Name='Euler');
        end

        function obj = rk4()
        % RK4 Klassisches Runge-Kutta-Verfahren (Ordnung 4).
            A = [0   0   0 0
                 1/2 0   0 0
                 0   1/2 0 0
                 0   0   1 0];
            b = [1/6 1/3 1/3 1/6];
            c = [0 1/2 1/2 1];
            obj = RungeKutta(A, b, c, Order=4, Name='RK4');
        end

        function obj = rkf45()
        % RKF45 Runge-Kutta-Fehlberg 4(5): propagiert mit Ordnung 5, Schätzer Ordnung 4.
            A = zeros(6);
            A(2,1)   = 1/4;
            A(3,1:2) = [3/32 9/32];
            A(4,1:3) = [1932/2197 -7200/2197 7296/2197];
            A(5,1:4) = [439/216 -8 3680/513 -845/4104];
            A(6,1:5) = [-8/27 2 -3544/2565 1859/4104 -11/40];
            c  = [0 1/4 3/8 12/13 1 1/2];
            b5 = [16/135 0 6656/12825 28561/56430 -9/50 2/55];
            b4 = [25/216 0 1408/2565 2197/4104 -1/5 0];
            obj = RungeKutta(A, b5, c, Order=5, bErr=b4, OrderErr=4, Name='RKF4(5)');
        end

        function obj = rkf78()
        % RKF78 Runge-Kutta-Fehlberg 7(8): propagiert mit Ordnung 8, Schätzer Ordnung 7.
            A = zeros(13);
            A(2,1)    = 2/27;
            A(3,1:2)  = [1/36 1/12];
            A(4,1:3)  = [1/24 0 1/8];
            A(5,1:4)  = [5/12 0 -25/16 25/16];
            A(6,1:5)  = [1/20 0 0 1/4 1/5];
            A(7,1:6)  = [-25/108 0 0 125/108 -65/27 125/54];
            A(8,1:7)  = [31/300 0 0 0 61/225 -2/9 13/900];
            A(9,1:8)  = [2 0 0 -53/6 704/45 -107/9 67/90 3];
            A(10,1:9) = [-91/108 0 0 23/108 -976/135 311/54 -19/60 17/6 -1/12];
            A(11,1:10) = [2383/4100 0 0 -341/164 4496/1025 -301/82 2133/4100 45/82 45/164 18/41];
            A(12,1:11) = [3/205 0 0 0 0 -6/41 -3/205 -3/41 3/41 6/41 0];
            A(13,1:12) = [-1777/4100 0 0 -341/164 4496/1025 -289/82 2193/4100 51/82 33/164 12/41 0 1];
            c  = [0 2/27 1/9 1/6 5/12 1/2 5/6 1/6 2/3 1/3 1 0 1];
            b8 = [0 0 0 0 0 34/105 9/35 9/35 9/280 9/280 0 41/840 41/840];
            b7 = [41/840 0 0 0 0 34/105 9/35 9/35 9/280 9/280 41/840 0 0];
            obj = RungeKutta(A, b8, c, Order=8, bErr=b7, OrderErr=7, Name='RKF7(8)');
        end

    end
end

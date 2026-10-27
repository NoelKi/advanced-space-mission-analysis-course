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

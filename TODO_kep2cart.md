# TODO: kep2cart in library.m korrigieren

## Syntax

- [ ] Rückgabewerte in eckige Klammern: `function [r_eci, v_eci] = kep2cart(kep, mu)`
- [ ] Drehmatrizen mit `;` zwischen den Zeilen schreiben (`[a b c; d e f; g h i]`), nicht `[[...],[...],[...]]` (das ergibt einen 1×9-Vektor)
- [ ] Perifokale Vektoren als Spaltenvektoren: `[cos(nu); sin(nu); 0]` statt `[cos(v), sin(v), 0]`
- [ ] `a(1-e^2)` → `a*(1-e^2)` (ohne `*` ist es ein Indexzugriff)
- [ ] `v` überall durch `nu` ersetzen
- [ ] `;` ans Ende der Zeilen 20–24 setzen, damit nichts im Command Window ausgegeben wird

## Physik / Mathematik

- [ ] Semi-latus rectum korrigieren: `p = a*(1 - e^2)` statt `a*(1-e)^2`
- [ ] `mu` definieren: als zweites Argument im `arguments`-Block, optional mit Standardwert
      `mu (1,1) double {mustBePositive} = 398600.4418   % km^3/s^2 (Erde)`
- [ ] Eine gemeinsame Drehmatrix für `r` **und** `v` verwenden:
      `Q = (R3_omega * R1_i * R3_RAAN)';   % perifokal -> inertial (transponiert!)`
- [ ] `orbitHeight` in `r` umbenennen (es ist der Abstand zum Erdmittelpunkt, nicht die Höhe)

## Skript

- [ ] `kep` mit echten Zahlen belegen, Winkel im Bogenmaß (`deg2rad(...)`)
- [ ] Funktion aufrufen: `[r, v] = kep2cart(kep);`
- [ ] Ausführbaren Code oben lassen, Funktionsdefinition(en) ans Dateiende

## Optional

- [ ] Position und Geschwindigkeit in zwei Funktionen trennen (`kep2pos`, `kep2vel`) mit gemeinsamer Hilfsfunktion `orbitBasics` für `p` und `Q`
- [ ] Wenn die Funktionen aus anderen Dateien aufgerufen werden sollen: jede in eine eigene `.m`-Datei mit gleichem Namen

## Tests

- [ ] Kreisbahn (e = 0): `norm(r) == a` und `norm(v) == sqrt(mu/a)`
- [ ] Alle Winkel = 0: `r == [a*(1-e); 0; 0]`
- [ ] Vis-Viva: `norm(v)^2/2 - mu/norm(r)` ergibt `-mu/(2*a)`
- [ ] Drehimpuls `h = cross(r, v)`: `acos(h(3)/norm(h))` ergibt wieder `i`
- [ ] Optional mit `keplerian2ijk` aus der Aerospace Toolbox vergleichen (falls vorhanden)

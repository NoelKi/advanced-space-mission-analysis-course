# MATLAB-Operatoren – Übersicht

Zu jedem Operator gibt es eine gleichwertige Funktion (z. B. `a + b` ist `plus(a, b)`).
Die steht in der Spalte **Funktion**; unter diesem Namen findest du ihn auch in der MATLAB-Hilfe (`doc plus`).

## 1. Arithmetische Operatoren

Faustregel: **mit Punkt = elementweise**, **ohne Punkt = Matrix-Rechnung** (lineare Algebra).
Bei `+` und `-` gibt es keinen Unterschied, deshalb haben sie keine Punkt-Version.

| Operator | Bedeutung | Funktion | Beispiel |
|---|---|---|---|
| `+` | Addition | `plus` | `[1 2] + [3 4]` → `[4 6]` |
| `-` | Subtraktion | `minus` | `[5 5] - [1 2]` → `[4 3]` |
| `+a` | Unäres Plus | `uplus` | `+3` → `3` |
| `-a` | Unäres Minus (Vorzeichen) | `uminus` | `-[1 2]` → `[-1 -2]` |
| `.*` | Elementweise Multiplikation | `times` | `[1 2] .* [3 4]` → `[3 8]` |
| `*` | Matrixmultiplikation | `mtimes` | `[1 2] * [3; 4]` → `11` |
| `./` | Elementweise Division (rechts) | `rdivide` | `[6 8] ./ [3 4]` → `[2 2]` |
| `.\` | Elementweise Division (links): `a .\ b = b ./ a` | `ldivide` | `[3 4] .\ [6 8]` → `[2 2]` |
| `/` | Matrix-Division rechts: `B / A ≈ B * inv(A)` | `mrdivide` | `x = b / A` löst `x*A = b` |
| `\` | Matrix-Division links: `A \ b ≈ inv(A) * b` | `mldivide` | `x = A \ b` löst `A*x = b` (Gleichungssystem!) |
| `.^` | Elementweise Potenz | `power` | `[1 2 3] .^ 2` → `[1 4 9]` |
| `^` | Matrixpotenz | `mpower` | `A^2` = `A*A` |
| `.'` | Transponieren | `transpose` | `[1 2 3].'` → Spaltenvektor |
| `'` | Konjugiert-komplex transponieren | `ctranspose` | bei reellen Zahlen gleich wie `.'` |

> **Typische Falle:** `r.^2` (jedes Element quadrieren) vs. `r^2` (Matrix mit sich selbst multiplizieren –
> bei einem 3x1-Vektor ein Fehler). Für Bahnradien also z. B. `vecnorm(r).^3`.

## 2. Vergleichsoperatoren

Ergebnis ist immer `logical` (`true`/`false` bzw. `1`/`0`), elementweise.

| Operator | Bedeutung | Funktion | Beispiel |
|---|---|---|---|
| `==` | gleich | `eq` | `[1 2 3] == 2` → `[0 1 0]` |
| `~=` | ungleich | `ne` | `[1 2 3] ~= 2` → `[1 0 1]` |
| `<` | kleiner | `lt` | `3 < 5` → `1` |
| `<=` | kleiner gleich | `le` | `5 <= 5` → `1` |
| `>` | größer | `gt` | `err > 1` |
| `>=` | größer gleich | `ge` | `tk >= tEnd` |

> **Achtung bei Gleitkommazahlen:** `0.1 + 0.2 == 0.3` ergibt `false`. Besser: `abs(a - b) < tol`.
> `isequal(A, B)` vergleicht ganze Arrays und liefert ein einziges `true`/`false`.

## 3. Logische Operatoren

| Operator | Bedeutung | Funktion | Beispiel |
|---|---|---|---|
| `~` | NICHT | `not` | `~isempty(x)` → „x ist nicht leer“ |
| `&` | UND, elementweise | `and` | `[1 0 1] & [1 1 0]` → `[1 0 0]` |
| `\|` | ODER, elementweise | `or` | `[1 0 1] \| [0 0 1]` → `[1 0 1]` |
| `&&` | UND mit Kurzschluss (nur Skalare) | – | `if ~isempty(b) && numel(b) == 3` |
| `\|\|` | ODER mit Kurzschluss (nur Skalare) | – | `if h < hMin \|\| err > 1` |
| – | Exklusives ODER | `xor` | `xor(true, true)` → `0` |
| – | mindestens ein Element wahr / alle wahr | `any`, `all` | `any(triu(A) ~= 0, 'all')` |

**Kurzschluss** (`&&`, `||`): Steht das Ergebnis nach der linken Seite schon fest, wird die rechte
gar nicht mehr ausgewertet. Deshalb ist `~isempty(x) && x(1) > 0` sicher – bei leerem `x` wird
`x(1)` nie aufgerufen. In `if`/`while`-Bedingungen nimmt man normalerweise `&&` und `||`.

## 4. Sonderzeichen und Syntax

| Zeichen | Bedeutung | Beispiel |
|---|---|---|
| `=` | Zuweisung | `a = 5` |
| `=` (im Aufruf) | Name-Wert-Argument (ab R2021a) | `newton(f, df, x0, Tol=1e-8)` |
| `:` | Bereich `start:schritt:ende` | `0:10:60` → `[0 10 20 ... 60]` |
| `:` (Index) | alle Elemente einer Dimension | `r(:, k)` → k-te Spalte; `A(:)` → alles als Spalte |
| `end` (Index) | letzter Index | `Y(:, end)`, `x(end+1) = 7` (anhängen) |
| `( )` | Indizierung, Funktionsaufruf, Klammerung | `A(2,3)`, `sin(x)`, `(a + b) * c` |
| `[ ]` | Arrays bauen, mehrere Rückgabewerte | `[1 2; 3 4]`, `[r, v] = kep2cart(kep)` |
| `[ ]` (leer) | leeres Array; löscht Elemente | `x(2) = []` |
| `{ }` | Cell-Arrays bauen bzw. Inhalt holen | `c = {1, 'abc'}`, `c{2}` → `'abc'` |
| `.` | Feld einer Struktur, Eigenschaft/Methode | `opts.Tol`, `OrbitKit.mu`, `rk.solve(...)` |
| `.( )` | dynamischer Feldname | `s.(name)` mit `name = 'Tol'` |
| `,` | trennt Argumente / Spalten / Befehle | `[1, 2, 3]`, `a = 1, b = 2` |
| `;` | neue Zeile im Array; unterdrückt Ausgabe | `[1; 2; 3]`, `a = 5;` |
| `...` | Zeile wird in der nächsten fortgesetzt | `error('id', ...` ⏎ `'Text')` |
| `%` | Kommentar | `% Umlaufzeit [s]` |
| `%%` | Abschnitt (Section) im Skript/Live Script | `%% Aufgabe 1b` |
| `%{ ... %}` | Blockkommentar (jeweils allein in der Zeile) | |
| `@` | Function Handle / anonyme Funktion | `f = @(t, y) -y`, `g = @sin` |
| `~` (Rückgabe) | Rückgabewert ignorieren | `[~, Y] = rk.solve(...)` |
| `~` (Eingabe) | Eingabeargument ignorieren | `function out = f(~, y)` |
| `'...'` | Zeichenvektor (`char`) | `'RKF4(5)'` |
| `"..."` | String (`string`) | `"Bahn"` |
| `!` | Befehl an das Betriebssystem | `!ls` |
| `?` | Metaklasse einer Klasse | `?OrbitKit` |

## 5. Rangfolge (Operator-Priorität)

Von **höchster** (wird zuerst ausgewertet) zu **niedrigster** Priorität. Gleiche Stufe: von links nach rechts.

1. Klammern `( )`
2. Transponieren und Potenz: `.'`, `'`, `.^`, `^`
3. Potenz mit unärem Vorzeichen/NICHT im Exponenten: `.^-`, `.^+`, `.^~`, `^-`, `^+`, `^~`
4. Unäres Plus/Minus und NICHT: `+a`, `-a`, `~a`
5. Multiplikation und Division: `.*`, `./`, `.\`, `*`, `/`, `\`
6. Addition und Subtraktion: `+`, `-`
7. Doppelpunkt `:`
8. Vergleiche: `<`, `<=`, `>`, `>=`, `==`, `~=`
9. Elementweises UND `&`
10. Elementweises ODER `|`
11. Kurzschluss-UND `&&`
12. Kurzschluss-ODER `||`

**Beispiele, die überraschen können:**

```matlab
-2^2        % = -4, weil ^ vor dem Vorzeichen kommt  ->  -(2^2)
2^-1        % = 0.5, Sonderregel aus Stufe 3
1:3 + 1     % = 1:4, weil + vor : kommt  ->  1:(3+1)
a < b < c   % NICHT wie in Mathe! Wird zu (a < b) < c, also 0 oder 1 < c
```

Im Zweifel: **Klammern setzen** – das ist auch für die Lesbarkeit besser.

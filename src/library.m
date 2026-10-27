% Library by Kieran Noel Mai for 
% (advanced) space mission analysis 

% Beispiel: a [km], e [-], Winkel in rad
kep = [7000, 0.01, deg2rad(98), deg2rad(30), deg2rad(45), deg2rad(10)];
[r_eci, v_eci] = SpaceMissionAnalysis.kep2cart(kep)

% Beispiel: Kepler-Gleichung M = E - e*sin(E) nach E lösen
e = 0.01;
M = deg2rad(10);
[E, iter] = SpaceMissionAnalysis.newton(@(E) E - e*sin(E) - M, @(E) 1 - e*cos(E), M)

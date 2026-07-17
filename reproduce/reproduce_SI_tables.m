function reproduce_SI_tables()
%REPRODUCE_SI_TABLES  Recompute SI Appendix Tables S3-S6 of Birch et al. (2023).
%
%   Runs every prediction behind SI Tables S3-S6 (and the main-text
%   Gale grain-size estimate) and prints the computed values next to the
%   values printed in the SI Appendix, so agreement can be checked line
%   by line. Uses the published Table S2 median constants throughout.
%
%   Conventions (verified to reproduce the published tables):
%     * Bedload columns: Eqs. 2/4/6, depth in the nH -> 0 form.
%     * Suspended columns: Eqs. 3/5/6 with D50 inverted from width and
%       slope (Eq. 5b). Again, I don't recall why this was chosen, but its
%       for you to decide now on what you think is best (I think in
%       retrospect, assuming a fixed D50 for sand is best). 
%     * Parenthetical ranges in Tables S3/S4 are 10th/90th-percentile
%       width envelopes at the central coefficients. Ranges in Tables
%       S5/S6 pair the (min width, min slope) and (max width, max slope)
%       corners.
%     * gamma = 0 everywhere except the Jezero suspended column
%       (gamma = 2); lambda = 0.35.
%
%   Known SI typos are flagged where I found them (again, sorry).

fprintf('================================================================\n');
fprintf(' Birch et al. (2023) PNAS -- SI table reproduction\n');
fprintf(' computed value   (printed value)\n');
fprintf('================================================================\n');

%% ---------------- Table S3: Peace Vallis Fan, Gale ---------------------
fprintf('\n--- Table S3: Peace Vallis Fan, Gale crater (Mars) ---\n');
fprintf('A = 730 km^2, V = 0.9 km^3, gamma = 0, B = 27 (11-43) m\n');
A = 730; V = 0.9; Bs = [27 11 43];

prn = { ...  % {slope, label, D50cm, Q, Qs, H, M, t}
  0.003, 'Distal fan  (S=0.003)', '4.8 (2.0-7.7)', '35 (3.7-112)', ...
      '1.3e-3 (1.4e-4-4e-3)', '1.5 (0.6-2.4)', '0.2 (0.02-0.6)', ''; ...
  0.01,  'Central fan (S=0.01) ', '26 (11-41)',    '60 (6.0-190)', ...
      '1e-2 (1e-3-0.04)',     '1.9 (0.8-3.0)', '0.3 (0.03-0.9)', ...
      '1.5e7 (4.6e6-1.4e8)'};

for i = 1:size(prn,1)
    o = predict_bedload(Bs, prn{i,1}, 'mars');
    m = formation_metrics(o.Q, o.Qs, 'area_km2', A, 'volume_km3', V, 'gamma', 0);
    fprintf('\n%s\n', prn{i,2});
    triple('  D50 (cm) ', o.D50*100, prn{i,3});
    triple('  Q (m3/s) ', o.Q,       prn{i,4});
    triple('  Qs (m3/s)', o.Qs,      prn{i,5});
    triple('  H (m)    ', o.H,       prn{i,6});
    triple('  M (mm/hr)', m.M,       prn{i,7});
    if ~isempty(prn{i,8})
        triple('  t (hours)', m.t,   prn{i,8});
    end
end

o = predict_bedload(Bs, 0.001, 'mars');
fprintf('\nMain text, S = 0.001:\n');
triple('  D50 (cm) ', o.D50*100, '1.0 +/- 0.6');

%% ---------------- Table S4: western Jezero delta -----------------------
fprintf('\n--- Table S4: western Jezero delta (Mars) ---\n');
fprintf('A = 12,700 km^2, V = 22.5 km^3, B = 45 (19-148) m\n');
A = 12700; V = 22.5; Bs = [45 19 148];

prn = { ...
  0.003, 'Bedload, low S  (S=0.003)', 0, '8 (3.4-26)',   '125 (15-2500)', ...
      '4.6e-3 (5e-4-9e-2)', '2.5 (1.1-8.3)', '3.6e-2 (4e-3-7e-1)', '8.8e8 (5e7-8e9)'; ...
  0.03,  'Bedload, high S (S=0.03) ', 0, '200 (85-660)', '357 (41-7000)', ...
      '2.7e-1 (3e-2-5.4)',  '3.8 (1.6-12.6)','1.0e-1 (1e-2-2.0)',  '1.5e7 (8e5-1e8)'};

for i = 1:size(prn,1)
    o = predict_bedload(Bs, prn{i,1}, 'mars');
    m = formation_metrics(o.Q, o.Qs, 'area_km2', A, 'volume_km3', V, 'gamma', prn{i,3});
    fprintf('\n%s   [gamma = %g]\n', prn{i,2}, prn{i,3});
    triple('  D50 (cm) ', o.D50*100, prn{i,4});
    triple('  Q (m3/s) ', o.Q,       prn{i,5});
    triple('  Qs (m3/s)', o.Qs,      prn{i,6});
    triple('  H (m)    ', o.H,       prn{i,7});
    triple('  M (mm/hr)', m.M,       prn{i,8});
    triple('  t (hours)', m.t,       prn{i,9});
end

o = predict_suspended(Bs, 0.003, 'mars', 'grain_size', 'invert');
m = formation_metrics(o.Q, o.Qs, 'area_km2', A, 'volume_km3', V, 'gamma', 2);
fprintf('\nSuspended load, low S (S=0.003)   [gamma = 2, D50 inverted]\n');
triple('  Q (m3/s) ', o.Q,  '230 (33-3200)');
triple('  Qs (m3/s)', o.Qs, '4e-3 (5e-4-3e-2)   [SI typo]');
triple('  H (m)    ', o.H,  '3.8 (2.0-9.1)');
triple('  M (mm/hr)', m.M,  '6.6e-2 (1e-2-9e-1)');
triple('  t (hours)', m.t,  '5.1e8 (5e7-3e9)');

%% ---------------- Table S5: Saraswati Flumen, Titan --------------------
fprintf('\n--- Table S5: Saraswati Flumen (Titan) ---\n');
fprintf('A = 10,000 km^2, V = 0.36 km^3, gamma = 0\n');
fprintf('corners: low = (B=175, S=4e-4), high = (B=700, S=2e-3)\n');
A = 10000; V = 0.36;
Bc = [175 700]; Sc = [4e-4 2e-3];      % corner-paired

cols = { ...
  'titan_cold', 'bed',  'cold, bedload  ', '1.6-61', '35-2300',  '6.5e-4-0.36', ...
      '1.0-5.6', '0.013-0.84', '1.8e5-1.8e8   [SI typo]'; ...
  'titan_warm', 'bed',  'warm, bedload  ', '1.5-56', '100-6900', '8.9e-4-0.49', ...
      '1.9-10.3','0.038-2.5',  '1.3e5-7.3e7'; ...
  'titan_cold', 'susp', 'cold, susp.load', 'N/A',    '680-4.5e4','1.2-9.4', ...
      '6.7-26',  '0.24-16',    '6.9e3-5.2e4'; ...
  'titan_warm', 'susp', 'warm, susp.load', 'N/A',    '730-4.9e4','0.43-3.3', ...
      '7.2-28',  '0.26-18',    '2.0e4-1.5e5'};

for i = 1:size(cols,1)
    if strcmp(cols{i,2}, 'bed')
        o = predict_bedload(Bc, Sc, cols{i,1});
    else
        o = predict_suspended(Bc, Sc, cols{i,1}, 'grain_size', 'invert');
    end
    m = formation_metrics(o.Q, o.Qs, 'area_km2', A, 'volume_km3', V, 'gamma', 0);
    fprintf('\nSaraswati [%s]\n', cols{i,3});
    if strcmp(cols{i,2}, 'bed')
        pair('  D50 (cm) ', o.D50*100, cols{i,4});
    end
    pair('  Q (m3/s) ', o.Q,  cols{i,5});
    pair('  Qs (m3/s)', o.Qs, cols{i,6});
    pair('  H (m)    ', o.H,  cols{i,7});
    pair('  M (mm/hr)', m.M,  cols{i,8});
    pair('  t (hours)', m.t,  cols{i,9});
end

%% ---------------- Table S6: Vid Flumina, Titan -------------------------
fprintf('\n--- Table S6: Vid Flumina (Titan) ---\n');
fprintf('A = 1,270 km^2  [main text prints 90,000 km^2, that aint right, this one is, sorry.');
fprintf('corners: low = (B=100, S=1.1e-3), high = (B=175, S=1.5e-3)\n');
A = 1270;
Bc = [100 175]; Sc = [1.1e-3 1.5e-3];

cols = { ...
  'titan_cold', 'bed',  'cold, bedload  ', '3.8-10',  '14-64',    '0.96e-3-6.7e-3', '0.72-1.3', '0.04-0.18'; ...
  'titan_warm', 'bed',  'warm, bedload  ', '3.5-9.4', '41-190',   '1.3e-3-9.2e-3',  '1.3-2.4',  '0.11-0.53'; ...
  'titan_cold', 'susp', 'cold, susp.load', 'N/A',     '405-1700', '0.26-0.68',      '5.5-8.9',  '1.1-4.8'; ...
  'titan_warm', 'susp', 'warm, susp.load', 'N/A',     '440-1900', '0.09-0.24',      '5.9-9.6',  '1.2-5.2'};

for i = 1:size(cols,1)
    if strcmp(cols{i,2}, 'bed')
        o = predict_bedload(Bc, Sc, cols{i,1});
    else
        o = predict_suspended(Bc, Sc, cols{i,1}, 'grain_size', 'invert');
    end
    m = formation_metrics(o.Q, o.Qs, 'area_km2', A);
    fprintf('\nVid Flumina [%s]\n', cols{i,3});
    if strcmp(cols{i,2}, 'bed')
        pair('  D50 (cm) ', o.D50*100, cols{i,4});
    end
    pair('  Q (m3/s) ', o.Q,  cols{i,5});
    pair('  Qs (m3/s)', o.Qs, cols{i,6});
    pair('  H (m)    ', o.H,  cols{i,7});
    pair('  M (mm/hr)', m.M,  cols{i,8});
end

end

% ------------------------------------------------------------------------
function triple(label, v, printed)
% Print bold value + (10th-90th width) range vs printed string.
fprintf('%s  %-32s (printed: %s)\n', label, fmt3(v), printed);
end

function pair(label, v, printed)
% Print low/high corner values vs printed range string.
fprintf('%s  %-32s (printed: %s)\n', label, fmt2(v), printed);
end

function s = fmt3(v)
s = sprintf('%s (%s-%s)', g3(v(1)), g3(v(2)), g3(v(3)));
end

function s = fmt2(v)
s = sprintf('%s-%s', g3(v(1)), g3(v(2)));
end

function s = g3(x)
s = sprintf('%.3g', x);
end

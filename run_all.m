%RUN_ALL  Set up paths, run a quick demo, and reproduce the SI tables.
%
%   From the repository root:
%       >> run_all
%
%   Runs in MATLAB (R2016b or later). No toolboxes are required.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));
addpath(fullfile(here, 'fitting'));
addpath(fullfile(here, 'reproduce'));

fprintf('birch2023-river-flows: quick demo\n');
fprintf('---------------------------------\n');
o = predict_bedload(27, 0.003, 'mars');
fprintf(['Gale distal fan (B = 27 m, S = 0.003, Mars, bedload):\n' ...
         '  D50 = %.2g cm, Q = %.0f m^3/s, H = %.1f m, Qs = %.2g m^3/s\n\n'], ...
        o.D50*100, o.Q, o.H, o.Qs);

reproduce_SI_tables();

fprintf(['\nOptional: refit the Table S2 constants from the compilation with\n' ...
         '  >> fits = fit_hydraulic_geometry();\n' ...
         '(bootstrap, ~1000 resamples; takes a minute or two).\n']);

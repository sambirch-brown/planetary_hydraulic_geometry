function c = published_coefficients()
%PUBLISHED_COEFFICIENTS  Fitted dimensionless hydraulic-geometry constants.
%
%   c = PUBLISHED_COEFFICIENTS() returns the empirical constants of
%   Table S2 in Birch et al. (2023), PNAS, doi:10.1073/pnas.2206837120,
%   as printed. These are the medians of the bootstrap distributions
%   (Materials & Methods, Sec. 5.2), with 16th/84th-percentile bounds
%   stored in .lo /.hi companions.
%
%   c.bedload    constants for bedload-dominated rivers (Eq. 1, a simple
%                power law), fit to gravel-bedded rivers with floodplains
%   c.suspended  constants for suspended load-dominated rivers (Eq. 1
%                with an extra Rep^m factor), fit to all sand-bedded rivers
%
%   c.D50_sand   median grain size of the suspended load-dominated rivers
%                in the terrestrial compilation [m]; used as the fixed
%                grain size when evaluating Rep-dependent terms for
%                planetary suspended-load predictions (main text, Eq. 3)
%
%   Regenerate with fitting/fit_hydraulic_geometry.m.

% ---- bedload-dominated (gravel with floodplains) -----------------------
b = struct();
b.alphaB = 5.2;   b.alphaB_lo = 4.2;   b.alphaB_hi = 6.2;
b.nB     = 0.06;  b.nB_lo     = 0.04;  b.nB_hi     = 0.08;
b.alphaH = 0.44;  b.alphaH_lo = 0.38;  b.alphaH_hi = 0.50;
b.nH     = -0.02; b.nH_lo     = -0.04; b.nH_hi     = 0.00;
b.alphaS = 0.09;  b.alphaS_lo = 0.06;  b.alphaS_hi = 0.12;
b.nS     = -0.33; b.nS_lo     = -0.36; b.nS_hi     = -0.30;

% ---- suspended load-dominated (all sand) -------------------------------
s = struct();
s.alphaB = 0.70;  s.alphaB_lo = 0.30;  s.alphaB_hi = 1.80;   % 0.70 -0.4/+1.1
s.nB     = 0.11;  s.nB_lo     = 0.08;  s.nB_hi     = 0.14;
s.mB     = 0.10;  s.mB_lo     = 0.01;  s.mB_hi     = 0.19;
s.alphaH = 3.7;   s.alphaH_lo = 2.1;   s.alphaH_hi = 6.0;    % 3.7 -1.6/+2.3
s.nH     = -0.06; s.nH_lo     = -0.08; s.nH_hi     = -0.04;
s.mH     = -0.11; s.mH_lo     = -0.18; s.mH_hi     = -0.04;
s.alphaS = 0.02;  s.alphaS_lo = 0.01;  s.alphaS_hi = 0.08;   % 0.02 -0.01/+0.06
s.nS     = -0.17; s.nS_lo     = -0.21; s.nS_hi     = -0.13;
s.mS     = -0.05; s.mS_lo     = -0.20; s.mS_hi     = 0.10;

c = struct('bedload', b, 'suspended', s, 'D50_sand', 3.45e-4);
end

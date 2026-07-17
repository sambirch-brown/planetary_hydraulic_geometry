function fits = fit_hydraulic_geometry(nboot)
%FIT_HYDRAULIC_GEOMETRY  Refit the dimensionless hydraulic-geometry constants.
%
%   fits = FIT_HYDRAULIC_GEOMETRY(nboot)
%
%   Repeats the bootstrap fitting of Birch et al. (2023), Materials &
%   Methods Sec. 5.2, against the terrestrial river compilation:
%
%     * Bedload-dominated rivers WITH floodplains (Type = 1,
%       WidthControl = 1): simple power laws in log-log space,
%           Bhat = alphaB Qhat^nB,  Hhat = alphaH Qhat^nH,
%           S    = alphaS Qhat^nS,
%       fit with a robust (bisquare / Tukey biweight) regression.
%
%     * All suspended load-dominated rivers (Type = 2): bivariate power
%       laws including the particle Reynolds number,
%           Bhat = alphaB Qhat^nB Rep^mB   (and similarly H, S),
%       fit by ordinary least squares in log space.
%
%   Each fit is bootstrapped nboot times (default 1000) by resampling
%   rivers with replacement; the reported constants are the medians of
%   the bootstrap distributions with 16th/84th-percentile bounds --
%   exactly the quantities printed in Table S2.
%
%   The random seed is fixed for repeatability of this script; bootstrap
%   medians will match Table S2 to within resampling noise.
%
%   Output: struct fits.bedload / fits.suspended with fields
%   alphaB, nB, (mB,) alphaH, nH, (mH,) alphaS, nS, (mS), each a struct
%   with .med / .lo / .hi, plus the raw bootstrap draws in .draws.

if nargin < 1, nboot = 1000; end
try, rng(42); catch, rand('state', 42); end

d = load_river_data();

% ---------------- bedload (gravel with floodplains) ---------------------
ig = d.Type == 1 & d.WidthControl == 1;
Qg = d.Q(ig); Bg = d.B(ig); Hg = d.H(ig); Sg = d.Slope(ig);
fprintf('Bedload-dominated, with floodplains: %d rivers\n', sum(ig));

fits.bedload = struct();
fits.bedload = fit_block(fits.bedload, 'B', log10(Qg), [], log10(Bg), nboot, true);
fits.bedload = fit_block(fits.bedload, 'H', log10(Qg), [], log10(Hg), nboot, true);
fits.bedload = fit_block(fits.bedload, 'S', log10(Qg), [], log10(Sg), nboot, true);

% ---------------- suspended (all sand) ----------------------------------
is = d.Type == 2 & ismember(d.WidthControl, [1 2 9]);
Qs = d.Q(is); Bs = d.B(is); Hs = d.H(is); Ss = d.Slope(is); Rs = d.Rp(is);
fprintf('Suspended load-dominated, all classes: %d rivers\n', sum(is));
Ds = d.D(is); Ds = Ds(isfinite(Ds));
fprintf('  median D50 of this class: %.3g m (fixed grain size of Eq. 3)\n', ...
        median(Ds));

fits.suspended = struct();
fits.suspended = fit_block(fits.suspended, 'B', log10(Qs), log10(Rs), log10(Bs), nboot, false);
fits.suspended = fit_block(fits.suspended, 'H', log10(Qs), log10(Rs), log10(Hs), nboot, false);
fits.suspended = fit_block(fits.suspended, 'S', log10(Qs), log10(Rs), log10(Ss), nboot, false);

% ---------------- report vs Table S2 ------------------------------------
c = published_coefficients();
fprintf('\n%-8s %-24s %s\n', 'param', 'this refit (med, 16-84)', 'Table S2');
rep('alphaB', fits.bedload.B.alpha,  c.bedload.alphaB);
rep('nB',     fits.bedload.B.n,      c.bedload.nB);
rep('alphaH', fits.bedload.H.alpha,  c.bedload.alphaH);
rep('nH',     fits.bedload.H.n,      c.bedload.nH);
rep('alphaS', fits.bedload.S.alpha,  c.bedload.alphaS);
rep('nS',     fits.bedload.S.n,      c.bedload.nS);
fprintf('  -- suspended --\n');
rep('alphaB', fits.suspended.B.alpha, c.suspended.alphaB);
rep('nB',     fits.suspended.B.n,     c.suspended.nB);
rep('mB',     fits.suspended.B.m,     c.suspended.mB);
rep('alphaH', fits.suspended.H.alpha, c.suspended.alphaH);
rep('nH',     fits.suspended.H.n,     c.suspended.nH);
rep('mH',     fits.suspended.H.m,     c.suspended.mH);
rep('alphaS', fits.suspended.S.alpha, c.suspended.alphaS);
rep('nS',     fits.suspended.S.n,     c.suspended.nS);
rep('mS',     fits.suspended.S.m,     c.suspended.mS);
end

% ------------------------------------------------------------------------
function fits = fit_block(fits, name, x1, x2, y, nboot, robust)
% Bootstrap fit of y = b0 + b1 x1 (+ b2 x2) in log space.
ok = isfinite(x1) & isfinite(y);
if ~isempty(x2), ok = ok & isfinite(x2); end
x1 = x1(ok); y = y(ok);
if ~isempty(x2), x2 = x2(ok); end
N = numel(y);

np = 2 + ~isempty(x2);
draws = zeros(nboot, np);
for i = 1:nboot
    idx = randi(N, N, 1);
    if isempty(x2), X = x1(idx); else, X = [x1(idx) x2(idx)]; end
    if robust
        b = robust_ls(X, y(idx));
    else
        b = [ones(N,1) X] \ y(idx);
    end
    draws(i,:) = b(:)';
end

q = @(v) struct('med', median(v), 'lo', prctile_(v,16), 'hi', prctile_(v,84));
out = struct('alpha', q(10.^draws(:,1)), 'n', q(draws(:,2)));
if np == 3, out.m = q(draws(:,3)); end
out.draws = draws;
fits.(name) = out;
end

function b = robust_ls(X, y)
% Iteratively reweighted least squares with a Tukey bisquare weight
% (tuning constant 4.685).
X = [ones(size(X,1),1) X];
b = X \ y;
for it = 1:50
    r = y - X*b;
    s = median(abs(r - median(r)))/0.6745;
    s = max(s, 1e-12);
    u = r/(4.685*s);
    w = (abs(u) < 1) .* (1 - u.^2).^2;
    sw = sqrt(w);
    bn = (X .* sw) \ (y .* sw);
    if norm(bn - b) < 1e-10*max(1, norm(b)), b = bn; return; end
    b = bn;
end
end

function p = prctile_(v, pct)
% Toolbox-free percentile (linear interpolation, as in prctile).
v = sort(v(:));
n = numel(v);
x = (0.5:n - 0.5)'/n*100;
p = interp1([0; x; 100], [v(1); v; v(end)], pct);
end

function rep(name, s, pub)
fprintf('%-8s %6.3f (%6.3f, %6.3f)   %g\n', name, s.med, s.lo, s.hi, pub);
end

# birch2023-river-flows

MATLAB code for **Birch et al. (2023), "Reconstructing river flows
remotely on Earth, Titan, and Mars," PNAS 120 (25), e2206837120**
([doi:10.1073/pnas.2206837120](https://doi.org/10.1073/pnas.2206837120)).

Given only a channel's bankfull width and bed slope, plus the fluid and
sediment properties of the planetary body, this code predicts the bed grain
size, bankfull discharge, flow depth, and sediment flux of an alluvial river,
using dimensionless hydraulic-geometry relations calibrated against a
compilation of terrestrial rivers. It reproduces every prediction in the
paper (SI Appendix, Tables S3–S6) and can refit the underlying empirical
constants (Table S2) from the included river compilation.

## Quick start

Runs in MATLAB (R2016b+) or GNU Octave (tested on 8.4). No toolboxes needed.

```matlab
>> run_all            % demo + full SI table reproduction
```

Or call the predictors directly:

```matlab
addpath src

% Bedload-dominated river (Gale crater distal fan, Mars):
out = predict_bedload(27, 0.003, 'mars');
%   out.D50 = 0.048 m,  out.Q = 35 m^3/s,  out.H = 1.5 m,  out.Qs = 1.3e-3 m^3/s

% Suspended load-dominated river (Titan, cold liquid), fixed sand D50:
out = predict_suspended(700, 2e-3, 'titan_cold');

% Runoff rate and deposit formation timescale (Eq. 9):
m = formation_metrics(out.Q, out.Qs, 'area_km2', 10000, ...
                      'volume_km3', 0.36, 'gamma', 0);
```

Planetary bodies: `'earth'`, `'mars'`, `'titan_cold'` (84 K, 0% C2H6),
`'titan_warm'` (91 K, 25% C2H6), or a custom struct with fields
`ps, rho, g, nu` in SI units (see `src/planet_properties.m`).

## Repository layout

```
run_all.m                        demo + reproduction driver
src/
  planet_properties.m            Table S1 physical constants
  published_coefficients.m       Table S2 fitted constants (medians, 16/84%)
  predict_bedload.m              Eqs. 2, 4a-c, 6  (gravel-bedded rivers)
  predict_suspended.m            Eqs. 3, 5a-c, 6  (sand-bedded rivers)
  formation_metrics.m            runoff M and formation time t (Eq. 9)
reproduce/
  reproduce_SI_tables.m          recomputes SI Tables S3-S6 line by line,
                                 printing computed vs published values
fitting/
  load_river_data.m              reads the terrestrial compilation
  fit_hydraulic_geometry.m       bootstrap refit of the Table S2 constants
data/
  FINAL_CLEANED_2023-2.txt       terrestrial river compilation (634 rivers)
```

## Brief method summary

The terrestrial calibration (Fig. 2, Table S2) fits dimensionless width,
depth, and slope against dimensionless discharge, separately for
bedload-dominated (gravel) rivers with floodplains (such that the channel 
width closure relation from Parker (1978) is most apt) and for suspended
load-dominated (sand) rivers; the sand fits carry an additional particle
Reynolds number (Rep) dependence. For other planetary bodies, the fitted
coefficients rescale with the submerged specific density R (Eqs. 2 and 3),
and the relations are inverted so that width and slope — the two quantities
measurable from orbit — predict D50, Q, H, and Qs (Eqs. 4–6).

Two conventions are worth knowing because they are what the published
numbers contain:

* **Depth uses the nH → 0 form of Eq. 4c** (H = αH Q^0.4 g^-0.2). The
  fitted nH is statistically indistinguishable from zero, consistent with
  the near-constant dimensionless depth found in Parker et al. (2007), so
  the simplified form is adopted. Every depth in Tables S3–S6 uses it.

* **The fixed sand grain size D50 = 3.45e-4 m** is the median D50 of the
  suspended load-dominated rivers in the compilation (verifiable with
  `fit_hydraulic_geometry`). It sets the Rep-dependent factors that become
  the planetary coefficients of Eq. 3, and it is the default grain size in
  `predict_suspended`.

## Reproducing the SI tables

`reproduce/reproduce_SI_tables.m` recomputes Tables S3–S6 and the
main-text Gale grain-size estimate, printing each computed value next to
the published one. Some notes: 

* **The Mars bedload entries (Tables S3, S4) should perfectly match what
  is in the paper** — bold central values and parenthetical ranges
  alike. The parenthetical ranges in Tables S3/S4 are 10th/90th-percentile
  *width* envelopes evaluated at the central coefficients; the ranges in
  Tables S5/S6 pair the (min width, min slope) and (max width, max slope)
  corners.

* **The Titan bedload entries (Tables S5, S6) reproduce for most quantities
  , with a few endpoints a tad off.** The published Saraswati column 
  evaluated some entries at the unrounded adjusted slopes
  (S = S_measured / (1.5 sin 22°) ≈ 3.6×10⁻⁴–2.1×10⁻³ rather than the
  rounded 4×10⁻⁴–2×10⁻³ quoted in the table), and the Titan predictions
  are a bit more sensitive than the Mars ones to the last printed digit of 
  the fitted constants (whoops). Re-running with the unrounded fit output  
  closes these gaps.

* **The suspended-load columns of Tables S4–S6 were generated with the
  grain size inverted from width and slope (Eq. 5b)**, not with the fixed
  median sand D50. Because the Eq. 5b inversion is sensitive to the last 
  printed digit of the fitted constantshis means a few points (especially 
  those with extreme discharges) are a bit off. Re-running the reproduction 
  with the original unrounded constants closes these gaps.

  **But note** Because the slope relation for sand rivers is relatively 
  insensitive to grain size, the inverted D50 itself is not meaningful for 
  suspended load-dominated rivers (the SI prints "N/A" for it), and the 
  recommended usage going forward is the fixed median sand grain size — 
  the default in `predict_suspended`. Pass `'grain_size','invert'` only 
  to reproduce the published tables (i.e., just use the default, and remember
  that my table values are a bit off, sorry).

Known typos in the published paper and SI (author-confirmed) are flagged
inline by the reproduction script.

## Refitting the empirical constants

```matlab
>> fits = fit_hydraulic_geometry();   % ~1000 bootstrap resamples
```

Bootstraps the log-space fits exactly as in Materials & Methods Sec. 5.2
(robust bisquare regression for the gravel relations, ordinary least
squares for the bivariate sand relations).

## Some typos in the paper/supplements I've found

I found some typos and internal inconsistencies in the published paper 
and SI Appendix while reformating all my messy code for this repository. 
None of them affect the paper's conclusions, and I documented each in the
code for you. If you find more, let me know (sambirch@brown.edu), I'm sure
there are a couple more :) 

1. **Table S5, cold-bedload formation time, upper bound.** Printed as
   1.8×10⁸ hours; the correct value is ~1.0×10⁸ hours (this code computes
   9.8×10⁷ h). The "1.8" was a copy error because I duplicated the 1.8×10⁵ h 
   lower bound. The corrected value is consistent with the main text, which
   states t > 1×10⁵–1×10⁸ h and ~2 Myr at the Titan storm intermittency.

2. **Vid Flumina drainage area.** Section 5.7.2 states A ~ 90,000 km²; the
   runoff rates in Table S6 (and the main-text range M = 0.04–5.2 mm/hr)
   were computed with A ≈ 1,270 km², which is the correct drainage area.
   The 90,000 km² figure in the text is an error. This repository uses
   A = 1,270 km².

3. **Table S4, suspended-load column, bold Qs.** The bold Qs = 4×10⁻³ m³/s
   is inconsistent with the bold t = 5.1×10⁸ h in the same column: with
   γ = 2 and λ = 0.35 (Eq. 9), that t implies Qs ≈ 2.7×10⁻³ m³/s, which is
   what this code computes. The parenthetical range values in the column
   are self-consistent.

4. **Rounding-level inconsistencies (because I'm lazy).**
   - The Eq. 2 exponent prints as 0.73, evaluated from the unrounded
     fitted Earth constants. Re-deriving it from the rounded Table S2
     medians gives 0.79, and likewise the prefactors 18/0.22/0.11/0.01
     become 17.7/0.20/0.118/0.0102 (~1–8% differences, negligible in reality
     given the other error sources). This repository uses Eq. 2 as printed, 
     which is what generated the published tables.
   - The Eq. 3 prefactor αB = 0.95 similarly reflects unrounded constants
     (rounded Table S2 values give 0.90).
   - The main-text Saraswati statement "Q = 40–50,000 m³/s" is the Table
     S5 envelope with loose rounding (the table minimum is 35 m³/s).

5. **Suspended-load grain-size convention.** The numerical values in the
   suspended-load columns of Tables S4–S6 were generated with D50 inverted
   from width and slope (Eq. 5b), not with the fixed median sand grain
   size (3.45×10⁻⁴ m) that defines the Eq. 3 coefficients. Since the
   inversion is unreliable for sand-bedded rivers (the SI prints "N/A" for
   suspended D50), the fixed grain size is the recommended convention
   going forward; the inversion is retained in `predict_suspended` under
   the option `'grain_size','invert'` solely to reproduce the published
   tables. See README.md, "Reproducing the SI tables."

## Citation

```bibtex
@article{birch2023river,
  author  = {Birch, S. P. D. and others},
  title   = {Reconstructing river flows remotely on Earth, Titan, and Mars},
  journal = {Proceedings of the National Academy of Sciences},
  volume  = {120},
  number  = {25},
  pages   = {e2206837120},
  year    = {2023},
  doi     = {10.1073/pnas.2206837120}
}
```

## License

MIT — see `LICENSE`.

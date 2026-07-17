function out = predict_suspended(B, S, planet, varargin)
%PREDICT_SUSPENDED  Bankfull predictions for a suspended load-dominated river.
%
%   out = PREDICT_SUSPENDED(B, S, planet)
%   out = PREDICT_SUSPENDED(B, S, planet, 'grain_size', gs)
%
%   Inputs
%       B       bankfull channel width [m]        (scalar or array)
%       S       channel-bed slope      [-]        (scalar or array)
%       planet  'earth' | 'mars' | 'titan_cold' | 'titan_warm', or a
%               custom struct (see planet_properties.m)
%
%   Option
%       'grain_size'  either a numeric bed grain size in meters, or the
%                     string 'invert' to solve for D50 from width and
%                     slope (Eq. 5b).
%                     Default: 3.45e-4 m, the median D50 of the
%                     suspended load-dominated rivers in the terrestrial
%                     compilation.
%
%   Output struct (elementwise over B, S)
%       out.D50    grain size used (fixed or inverted)   [m]
%       out.Q      bankfull discharge                    [m^3/s]
%       out.H      bankfull depth                        [m]
%       out.Qs     total sediment flux                   [m^3/s]
%       out.Rep    particle Reynolds number at out.D50
%       out.planet planetary properties used
%
%   Method (Birch et al. 2023, PNAS, doi:10.1073/pnas.2206837120):
%     The suspended load-dominated fits (Table S2) include an explicit
%     particle-Reynolds-number dependence, B~ = alphaB Qhat^nB Rep^mB
%     (and similarly for H~ and S). Writing the dimensional relations
%     out and letting Rep carry the (R, g, nu, D50) dependence gives
%
%       Q  = [ (B/alphaB) R^(-mB/2) g^(0.2+(nB-mB)/2)
%              D50^(2.5 nB - 1.5 mB) nu^mB ]^(1/(nB+0.4))          (Eq. 5a)
%       H  = alphaH Q^(nH+0.4) R^(mH/2) g^(-0.2+(mH-nH)/2)
%              D50^(1.5 mH - 2.5 nH) nu^(-mH)                       (Eq. 5c)
%       Qs = alphaY(Rep) Q^nY g^((1-nY)/2) D50^(2.5(1-nY))          (Eq. 6)
%
%     with nY = 0.4 - (nB + 3 nH + nS) + 2.5(0.4 + nH + nS) and
%     alphaY(Rep) built from the fitted constants and the Parker-style
%     prefactor (SI Appendix). Absorbing the Rep factors at fixed
%     D50 = 3.45e-4 m is what produces the planetary coefficients quoted
%     in main-text Eq. 3.
%
%     GRAIN SIZE. The paper notes that the slope relation for suspended
%     load-dominated rivers is relatively insensitive to grain size, so
%     inverting D50 from width and slope (Eq. 5b) is unreliable for
%     these rivers (it is why SI Tables S4-S6 print "N/A" for suspended
%     D50.) The recommended (default) usage therefore fixes D50 at the 
%     median sand grain size of the terrestrial dataset. With a fixed D50 
%     the predictions depend on width but not on slope.
%
%     Note, however, that the *numerical values* printed in the
%     suspended-load columns of SI Tables S4-S6 were generated with the
%     Eq. 5b inversion active (their ranges span slope as well as width,
%     which a fixed-D50 calculation cannot do). Use
%     'grain_size','invert' to reproduce the printed tables; see
%     README.md, "Reproducing the SI tables", for the comparison. It was
%     long enough ago now that there was a reason I did it this way, I just
%     don't recall why. Nevertheless, the default is the safest way to go
%     if you are doing your own work and not reproducing mine!
%
%   Example (default, fixed grain size):
%       out = predict_suspended(700, 2e-3, 'titan_cold');

gs = 3.45e-4;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'grain_size', gs = varargin{k+1};
        otherwise, error('Unknown option ''%s''.', varargin{k});
    end
end

p = planet_properties(planet);
c = published_coefficients();
s = c.suspended;

ab = s.alphaB;  nb = s.nB;  mb = s.mB;
ah = s.alphaH;  nh = s.nH;  mh = s.mH;
as = s.alphaS;  ns = s.nS;  ms = s.mS;
aa = 0.05;                       % entrainment prefactor (SI Appendix)

R  = p.R;  g = p.g;  nu = p.nu;

% ---- grain size: fixed (default) or inverted from width & slope --------
if ischar(gs) && strcmpi(gs, 'invert')
    D = ( (ab./B) .* (S./as).^((nb + 0.4)/ns) .* ...
          ((R*g/nu^2).^((-0.5*nb*(ms/ns)) + (0.5*mb) - (0.2*(ms/ns)))) ...
        ).^( 1/((1.5*nb*(ms/ns)) + (0.6*(ms/ns)) - (1.5*mb) - 1) );   % Eq. 5b
else
    D = gs .* ones(size(B + S));   % broadcast to common shape
end

% ---- predictions (Eqs. 5a, 5c, 6) --------------------------------------
Q  = ( (B./ab) .* R^(-0.5*mb) .* g^(0.2 + 0.5*(nb - mb)) .* ...
       D.^(2.5*nb - 1.5*mb) .* nu^mb ).^(1/(nb + 0.4));               % Eq. 5a
H  = ah .* Q.^(nh + 0.4) .* R^(0.5*mh) .* g^(-0.2 + 0.5*(mh - nh)) .* ...
     D.^(1.5*mh - 2.5*nh) .* nu^(-mh);                                % Eq. 5c

Rep = sqrt(R*g)/nu .* D.^1.5;
nY  = 0.4 - (nb + 3*nh + ns) + 2.5*(0.4 + nh + ns);
aY  = (aa/R^2) * (as^1.5/(ah^0.5 * ab)) .* Rep.^(1.5*ms - 0.5*mh - mb);
Qs  = aY .* Q.^nY .* g^(0.5*(1 - nY)) .* D.^(2.5*(1 - nY));           % Eq. 6

out = struct('D50', D, 'Q', Q, 'H', H, 'Qs', Qs, 'Rep', Rep, ...
             'planet', p);
end

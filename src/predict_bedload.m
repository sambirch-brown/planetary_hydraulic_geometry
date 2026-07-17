function out = predict_bedload(B, S, planet)
%PREDICT_BEDLOAD  Bankfull predictions for a bedload-dominated river.
%
%   out = PREDICT_BEDLOAD(B, S, planet)
%
%   Inputs
%       B       bankfull channel width [m]        (scalar or array)
%       S       channel-bed slope      [-]        (scalar or array)
%       planet  'earth' | 'mars' | 'titan_cold' | 'titan_warm', or a
%               custom struct (see planet_properties.m)
%
%   Output struct (elementwise over B, S)
%       out.D50    median bed grain size            [m]
%       out.Q      bankfull discharge               [m^3/s]
%       out.H      bankfull depth                   [m]
%       out.Qs     bedload sediment flux            [m^3/s]
%       out.alpha  planetary coefficients (alphaB, alphaH, alphaS, alphaY)
%       out.planet planetary properties used
%
%   Method (Birch et al. 2023, PNAS, doi:10.1073/pnas.2206837120):
%     The Earth-calibrated dimensionless relations (Table S2, bedload
%     column) are rescaled to the target body's submerged specific
%     density R = (ps - rho)/rho through the threshold-channel
%     derivation of the SI Appendix, which yields main-text Eq. 2:
%
%         alphaB = 18 / ((1+R) sqrt(R))
%         alphaH = 0.22 (1+R)^0.73
%         alphaS = 0.11 R (1+R)^(-0.73)
%         alphaY = 0.01 / (1+R)
%
%     The prefactors and the 0.73 exponent were evaluated using the
%     unrounded fitted Earth constants; they are used here as printed
%     because they are what generated every value in SI Tables S3-S6. 
%     Re-deriving them from the rounded Table S2 medians gives
%     17.7, 0.20, 0.118, 0.0102, and exponent 0.79 -- differences of a
%     few percent that don't really matter in reality. 
%
%     Predictions then follow Eqs. 4b (D50), 4a (Q), 4c (H), and 6 (Qs):
%         D50 = (B/alphaB) (alphaS/S)^((nB+0.4)/nS)
%         Q   = [ g^(0.2+nB/2) B D50^(2.5 nB) / alphaB ]^(1/(nB+0.4))
%         H   = alphaH Q^0.4 g^(-0.2)
%         Qs  = alphaY Q^nY g^((1-nY)/2) D50^(2.5(1-nY)),  nY = 1+nB+1.5 nS
%
%     The depth uses the nH -> 0 form of Eq. 4c. The fitted nH = -0.02
%     is statistically indistinguishable from zero (Table S2), matching
%     the near-constant dimensionless depth found by Parker et al.
%     (2007), so nH = 0 is adopted for simplicity and consistency with
%     that work. Every depth in SI Tables S3-S6 uses this form.
%
%   Example (Gale distal fan, Table S3):
%       out = predict_bedload(27, 0.003, 'mars');
%       % out.D50 = 0.048 m, out.Q = 35 m^3/s, out.H = 1.5 m

p = planet_properties(planet);
c = published_coefficients();
b = c.bedload;

nb = b.nB;
ns = b.nS;

% ---- planetary coefficients: main-text Eq. 2, as printed ---------------
R  = p.R;
ex = 0.73;
aB = 18   / ((1 + R)*sqrt(R));
aH = 0.22 * (1 + R)^ex;
aS = 0.11 * R*(1 + R)^(-ex);
aY = 0.01 / (1 + R);

% ---- predictions (Eqs. 4 and 6) ----------------------------------------
g  = p.g;
D  = (B./aB) .* (aS./S).^((nb + 0.4)/ns);                          % Eq. 4b
Q  = (g.^(0.2 + 0.5*nb) .* B .* D.^(2.5*nb) ./ aB).^(1/(nb + 0.4));% Eq. 4a
H  = aH .* Q.^0.4 .* g.^(-0.2);                                    % Eq. 4c
nY = 1 + nb + 1.5*ns;
Qs = aY .* Q.^nY .* g.^((1 - nY)/2) .* D.^(2.5*(1 - nY));          % Eq. 6

out = struct('D50', D, 'Q', Q, 'H', H, 'Qs', Qs, ...
             'alpha', struct('alphaB', aB, 'alphaH', aH, ...
                             'alphaS', aS, 'alphaY', aY, 'exponent', ex), ...
             'planet', p);
end

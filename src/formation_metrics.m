function out = formation_metrics(Q, Qs, varargin)
%FORMATION_METRICS  Runoff rate and deposit formation timescale.
%
%   out = FORMATION_METRICS(Q, Qs, 'area_km2', A, 'volume_km3', V, ...)
%
%   Inputs
%       Q    bankfull discharge      [m^3/s]  (scalar or array)
%       Qs   sediment flux           [m^3/s]  (scalar or array)
%
%   Options
%       'area_km2'        drainage area A [km^2]        -> enables out.M
%       'volume_km3'      deposit volume V [km^3]       -> enables out.t
%       'gamma'           washload fraction gamma        (default 0)
%       'lambda'          deposit porosity lambda        (default 0.35)
%       'hours_per_year'  flow hours per Earth year      -> enables out.t_yr
%                         (Mars, I = 0.1:  36.25*24 = 870;
%                          Titan, 2 storms per 2 Titan-day summer:
%                          (32*2*24)/28 = 54.86)
%
%   Output struct
%       out.M     runoff rate, Q/A                         [mm/hr]
%       out.t     active bankfull flow time to build the deposit [hours]
%                 t = (1-lambda)/(1+gamma) * V/Qs          (Eq. 9)
%       out.t_yr  t converted to Earth years at the given intermittency
%
%   Follows Birch et al. (2023), PNAS, doi:10.1073/pnas.2206837120.

A = []; V = []; gam = 0; lam = 0.35; hpy = [];
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'area_km2',       A   = varargin{k+1};
        case 'volume_km3',     V   = varargin{k+1};
        case 'gamma',          gam = varargin{k+1};
        case 'lambda',         lam = varargin{k+1};
        case 'hours_per_year', hpy = varargin{k+1};
        otherwise, error('Unknown option ''%s''.', varargin{k});
    end
end

out = struct();
if ~isempty(A)
    out.M = Q ./ (A*1e6) * 3600 * 1000;          % [mm/hr]
end
if ~isempty(V)
    out.t = (1 - lam)/(1 + gam) .* (V*1e9 ./ Qs) / 3600;   % [hours]
    if ~isempty(hpy)
        out.t_yr = out.t / hpy;                  % [Earth years]
    end
end
end

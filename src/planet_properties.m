function p = planet_properties(body)
%PLANET_PROPERTIES  Physical constants for Earth, Mars, and Titan (Table S1).
%
%   p = PLANET_PROPERTIES(body) returns a struct with fields:
%       ps  - sediment density        [kg/m^3]
%       rho - fluid density           [kg/m^3]
%       g   - gravity                 [m/s^2]
%       nu  - kinematic viscosity     [m^2/s]
%       R   - submerged specific density, (ps - rho)/rho
%
%   body is one of:
%       'earth'       quartz sediment, liquid water
%       'mars'        basaltic sediment, liquid water
%       'titan_cold'  water-ice sediment, 84 K liquid, 0% C2H6
%       'titan_warm'  water-ice sediment, 91 K liquid, 25% C2H6
%
%   Alternatively, pass a struct with fields ps, rho, g, nu (SI units)
%   to use custom properties; R is filled in automatically.
%
%   Values follow Table S1 of Birch et al. (2023), PNAS,
%   doi:10.1073/pnas.2206837120.

if isstruct(body)
    p = body;
    if ~all(isfield(p, {'ps','rho','g','nu'}))
        error('Custom planet struct needs fields ps, rho, g, nu (SI units).');
    end
    p.R = (p.ps - p.rho)/p.rho;
    if ~isfield(p, 'name'), p.name = 'custom'; end
    return
end

switch lower(body)
    case 'earth'
        p = struct('ps', 2650, 'rho', 1000, 'g', 9.81,  'nu', 0.01e-4);
    case 'mars'
        p = struct('ps', 2900, 'rho', 1000, 'g', 3.711, 'nu', 0.01e-4);
    case 'titan_cold'
        p = struct('ps', 950,  'rho', 670,  'g', 1.35,  'nu', 0.003e-4);
    case 'titan_warm'
        p = struct('ps', 950,  'rho', 540,  'g', 1.35,  'nu', 0.006e-4);
    otherwise
        error('Unknown body ''%s''. Use earth, mars, titan_cold, or titan_warm.', body);
end
p.R = (p.ps - p.rho)/p.rho;
p.name = lower(body);
end

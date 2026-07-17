function d = load_river_data(fname)
%LOAD_RIVER_DATA  Reads the terrestrial river compilation.
%
%   d = LOAD_RIVER_DATA(fname) reads the tab-separated compilation
%   (data/FINAL_CLEANED_2023-2.txt) and returns a struct with one field
%   per column:
%
%     B, H, Q      dimensionless width, depth, discharge (Bhat, Hhat, Qhat)
%     Rp           particle Reynolds number
%     Slope        channel-bed slope
%     Width, Depth, Discharge   dimensional values [m], [m], [m^3/s]
%     D            bed grain size D50 [m]
%     SedQ, SedQ_err            sediment flux and its uncertainty [m^3/s]
%     Latitude, Longitude
%     Type         1 = bedload-dominated (gravel), 2 = suspended (sand)
%     WidthControl 1 = free (floodplains), 2 = fixed/confined, other =
%                  additional classes (3, 9, ...) as in the compilation
%
%   Missing entries are NaN.

if nargin < 1
    here  = fileparts(mfilename('fullpath'));
    fname = fullfile(here, '..', 'data', 'FINAL_CLEANED_2023-2.txt');
end

fid = fopen(fname, 'r');
if fid < 0, error('Cannot open %s', fname); end
hdr   = fgetl(fid);
names = strsplit(strtrim(hdr), sprintf('\t'));
raw   = textscan(fid, repmat('%f', 1, numel(names)), 'Delimiter', sprintf('\t'));
fclose(fid);

d = struct();
for i = 1:numel(names)
    d.(names{i}) = raw{i};
end
d.n = numel(d.(names{1}));
end

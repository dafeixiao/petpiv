function frame = read_B_snapshot(fname, iflargefile)
% read_B_snapshot
% -------------------------------------------------------------------------
% Reads ONE Tecplot ASCII POINT snapshot (one ZONE) from a B*.dat file.
% Works for files WITH or WITHOUT uncertainty columns.
%
% - Detects VARIABLES whether QUOTED or UNQUOTED.
% - If uncertainty columns missing -> fill with UNC_DEFAULT = 1e-6.
% - Outputs consistent fields for downstream CCM/SVD solvers.
% -------------------------------------------------------------------------

UNC_DEFAULT = 1e-6;

if nargin < 2, iflargefile = 0; end
if iflargefile ~= 0
    error(['read_B_snapshot: iflargefile=1 not supported in this auto-detect version.\n' ...
           'Keep the VARIABLES line (iflargefile=0), or ask for fixed-column reader.']);
end

fid = fopen(fname,'r');
if fid < 0, error('Cannot open file: %s', fname); end
c = onCleanup(@() fclose(fid));

% --- TITLE
tline = fgetl(fid);
if ~ischar(tline), error('Empty file: %s', fname); end 

% --- VARIABLES line
vline = fgetl(fid);
if ~ischar(vline) || ~contains(vline, 'VARIABLES', 'IgnoreCase', true)
    error('Missing or invalid VARIABLES line in: %s', fname);
end

% Parse variable names robustly (quoted OR unquoted)
varnames = parse_tecplot_variables_line(vline);
M = numel(varnames);
if M < 3
    error('Parsed too few VARIABLES (%d) in %s. VARIABLES line: %s', M, fname, vline);
end

% --- Find I= ... (point count)
N = NaN;
while true
    line = fgetl(fid);
    if ~ischar(line)
        error('Could not find I=... line in %s', fname);
    end
    tI = regexp(line, 'I\s*=\s*(\d+)\s*,', 'tokens', 'once');
    if ~isempty(tI)
        N = str2double(tI{1});
        break;
    end
end
if ~isfinite(N) || N <= 0, error('Invalid point count parsed in %s', fname); end

% --- Skip header lines until numeric starts
pos = ftell(fid);
ln  = fgetl(fid);
while ischar(ln) && isempty(regexp(strtrim(ln), '^[\+\-]?\d', 'once'))
    pos = ftell(fid);
    ln  = fgetl(fid);
end
fseek(fid, pos, 'bof');

% --- Read numeric block
data = fscanf(fid, '%f', [M, N]).';
if size(data,1) ~= N || size(data,2) ~= M
    error('Numeric read mismatch in %s (got %dx%d, expected %dx%d)', ...
        fname, size(data,1), size(data,2), N, M);
end

getcol = @(names, defaultVal) get_col_by_name(data, varnames, names, defaultVal, N);

% -------------------------------------------------------------------------
% REQUIRED POSITION (your file uses x[mm], y[mm], z[mm])
% -------------------------------------------------------------------------
frame.X = get_required(getcol(["x","X","x[mm]","X[mm]","x (mm)","X (mm)","x_mm","X_mm"], NaN), ...
    "X", fname, varnames);

frame.Y = get_required(getcol(["y","Y","y[mm]","Y[mm]","y (mm)","Y (mm)","y_mm","Y_mm"], NaN), ...
    "Y", fname, varnames);

frame.Z = get_required(getcol(["z","Z","z[mm]","Z[mm]","z (mm)","Z (mm)","z_mm","Z_mm"], NaN), ...
    "Z", fname, varnames);

% -------------------------------------------------------------------------
% VELOCITY (your file uses Vx[m/s], Vy[m/s], Vz[m/s])
% -------------------------------------------------------------------------
frame.U = getcol(["U","u","Vx","vx","Vx[m/s]","vx[m/s]","Vx (m/s)","VelX","VelocityX"], 0);
frame.V = getcol(["V","v","Vy","vy","Vy[m/s]","vy[m/s]","Vy (m/s)","VelY","VelocityY"], 0);
frame.W = getcol(["W","w","Vz","vz","Vz[m/s]","vz[m/s]","Vz (m/s)","VelZ","VelocityZ"], 0);

% -------------------------------------------------------------------------
% ACCELERATION (your file uses Ax[m/s²], Ay[m/s²], Az[m/s²])
% Note: some exports use '^2' instead of the unicode ², so we include both.
% -------------------------------------------------------------------------
frame.ax = getcol(["ax","Ax","Ax[m/s²]","Ax[m/s^2]","Ax (m/s^2)","accx","AccX","AccelerationX"], 0);
frame.ay = getcol(["ay","Ay","Ay[m/s²]","Ay[m/s^2]","Ay (m/s^2)","accy","AccY","AccelerationY"], 0);
frame.az = getcol(["az","Az","Az[m/s²]","Az[m/s^2]","Az (m/s^2)","accz","AccZ","AccelerationZ"], 0);

% -------------------------------------------------------------------------
% TRACK ID
% -------------------------------------------------------------------------
frame.ID = getcol(["trackID","TrackID","ID","Id","id","ParticleID"], 0);

% -------------------------------------------------------------------------
% UNCERTAINTIES (if absent -> UNC_DEFAULT)
% Velocity uncertainty
frame.UNCU  = getcol(["UncertaintyVx","UncVx","UNCU","uncu","UNC_Vx","unc_Vx"], UNC_DEFAULT);
frame.UNCV  = getcol(["UncertaintyVy","UncVy","UNCV","uncv","UNC_Vy","unc_Vy"], UNC_DEFAULT);
frame.UNCW  = getcol(["UncertaintyVz","UncVz","UNCW","uncw","UNC_Vz","unc_Vz"], UNC_DEFAULT);

% Acceleration uncertainty
frame.UNCAx = getcol(["UncertaintyAx","UncAx","UNCAx","uncax","UNC_Ax","unc_Ax"], UNC_DEFAULT);
frame.UNCAy = getcol(["UncertaintyAy","UncAy","UNCAy","uncay","UNC_Ay","unc_Ay"], UNC_DEFAULT);
frame.UNCAz = getcol(["UncertaintyAz","UncAz","UNCAz","uncaz","UNC_Az","unc_Az"], UNC_DEFAULT);

% (Optional position uncertainty – keep consistent outputs)
frame.UNCX  = getcol(["UncertaintyX","UncX","UNCX","uncx"], UNC_DEFAULT);
frame.UNCY  = getcol(["UncertaintyY","UncY","UNCY","uncy"], UNC_DEFAULT);
frame.UNCZ  = getcol(["UncertaintyZ","UncZ","UNCZ","uncz"], UNC_DEFAULT);

end

% =========================================================================
% Helpers
% =========================================================================
function varnames = parse_tecplot_variables_line(vline)
% Handles:
%   VARIABLES = "x","y","z",...
% and also:
%   VARIABLES = x[mm], y[mm], z[mm], ...

% First try quoted parsing
tok = regexp(vline, '"([^"]+)"', 'tokens');
if ~isempty(tok)
    varnames = string([tok{:}]);
    varnames = strip(varnames);
    return;
end

% Fallback: unquoted list parsing
% Remove leading 'VARIABLES' and everything up to '='
eqpos = strfind(vline, '=');
if isempty(eqpos)
    % If weird formatting, just split after 'VARIABLES'
    tmp = regexprep(vline, '(?i)VARIABLES', '');
else
    tmp = vline(eqpos(1)+1:end);
end

parts = split(string(tmp), ',');
parts = strip(parts);

% Remove any accidental surrounding quotes
parts = replace(parts, '"', '');
parts = replace(parts, '''', '');

% Remove empty entries
parts(parts == "") = [];

varnames = parts;
end

function col = get_col_by_name(data, varnames, candidateNames, defaultVal, N)
idx = [];
for k = 1:numel(candidateNames)
    j = find(strcmpi(varnames, candidateNames(k)), 1, 'first');
    if ~isempty(j)
        idx = j;
        break;
    end
end
if isempty(idx)
    col = defaultVal * ones(N,1);
else
    col = data(:,idx);
end
end

function out = get_required(col, name, fname, varnames)
if all(isnan(col))
    error('Required position variable "%s" not found in %s. VARIABLES were: %s', ...
        name, fname, strjoin(string(varnames), ', '));
end
out = col;
end

%% CCM2 + SVD batch 

clear; clc;

%% ---------------- User parameters ----------------
path = '.\dat_output\';
iflargefile = 0;
% Choose processing range (center frames)
Tstart = 2;          % first frame to process. Make sure there is atleast one frame before it
Tend   = 2;         % last frame to process. Make sure there is atleast one frame after it

%% Output
dir_save_ccm = fullfile(path, 'CCM output');  %  Save directory 
dir_save_svd = fullfile(path, 'SVD output');  %  Save directory
ifsave_svd   = 1;                             %  flag to save svd or not

fileName_prefix_SVD = '';
fileName_prefix_CCM = '';

%% Averaging cell size (half-size)
Lx = 2; Ly = 2; Lz = 2;
Region = [Lx, Ly, Lz];

%% Minimum number of unstructured vectors
Nm_v = 10;
Nm_a = 10;

%% Grid parameters
% delta      = 2;      % Grid resolution in mm
delta      = 3;      % Grid resolution in mm
scale      = 1;
nu         = 1e-6;   % Kinematic viscosity of water
scale_unit = 1e-3;   % mm to m

XLimits=[-20 12];    %  X Grid limits
YLimits=[-24 24];    %  Y Grid limits

ZLimits=[0 20];     %  Z Grid limits

Xg = XLimits(1)+0.5*delta : 0.5*delta : XLimits(2)-0.5*delta;
Yg = YLimits(1)+0.5*delta : 0.5*delta : YLimits(2)-0.5*delta;
Zg = ZLimits(1)+0.5*delta : 0.5*delta : ZLimits(2)-0.5*delta;

Nz = numel(Zg); cutoff = ones(Nz,1);
if Nz>=7, cutoff([1:2 end-1:end]) = 2.5;
elseif Nz>=5, cutoff([1 end]) = 2.5; end
cutoff_smooth = smoothdata(cutoff,"gaussian",3)*0.01;

N1 = numel(Xg); N2 = numel(Yg); N3 = numel(Zg);
totalN = N1*N2*N3;

plotline1  = 'TITLE = "INITIAL"';
plotline2v = 'VARIABLES = "x","y","z","U","V","W","dudx","dudy","dudz","dvdx","dvdy","dvdz","dwdx","dwdy","dwdz","isValid"';
plotline2a = 'VARIABLES = "x","y","z","ax","ay","az","isValid"';
plotline3  = ['ZONE T="zeros" I=', num2str(N1,'%d'), ', J=', num2str(N2,'%d'), ', K=', num2str(N3,'%d'), ', F=POINT'];

ensure_dir(dir_save_ccm);
if ifsave_svd, ensure_dir(dir_save_svd); end

%% ---------------- List and sort B files ----------------
if size(dir(fullfile(path, 'B*.dat')), 1) ~= 0
    B = dir(fullfile(path, 'B*.dat'));
    if isempty(B), error('No B*.dat files in %s', path); end
    [~,ix] = sort({B.name}); B = B(ix);
    N_files = numel(B);
    
    % Map filename numeric ID: B000007.dat -> 7
    Bnum = nan(N_files,1);
    for i=1:N_files
        tok = regexp(B(i).name,'B(\d+)\.dat','tokens','once');
        if ~isempty(tok), Bnum(i) = str2double(tok{1}); end
    end
else
    B = dir(fullfile(path, 'TecPlot*.dat'));
    if isempty(B), error('No TecPlot*.dat files in %s', path); end
    [~,ix] = sort({B.name}); B = B(ix);
    N_files = numel(B);
    
    Bnum = nan(N_files,1);
    for i=1:N_files
        tok = regexp(B(i).name,'TecPlot(\d+)\.dat','tokens','once');
        if ~isempty(tok), Bnum(i) = str2double(tok{1}); end
    end
end

%% ---------------- Interpret Tstart/Tend ----------------
% If user typed values that match filename numbers, use those.
% Else treat them as file indices.
iStart = find(Bnum==Tstart, 1, 'first');
iEnd   = find(Bnum==Tend,   1, 'first');

if isempty(iStart) || isempty(iEnd)
    % fallback: treat as indices
    iStart = Tstart;
    iEnd   = Tend;
end

% Safety: need neighbors for 3-frame window
if iStart < 2
    error('Tstart too early. Need previous frame. Choose Tstart>=2 (center index).');
end
if iEnd > (N_files-1)
    error('Tend too late. Need next frame. Choose Tend<=N_files-1 (center index).');
end
if iEnd < iStart
    error('Tend must be >= Tstart.');
end

fprintf('Total files: %d\n', N_files);
fprintf('Processing CENTER range: iStart=%d (%s)  to  iEnd=%d (%s)\n', ...
    iStart, B(iStart).name, iEnd, B(iEnd).name);

%% ---------------- Initialize window at iStart ----------------
% Need (iStart-1, iStart, iStart+1)
frame_last = read_B_snapshot(fullfile(path, B(iStart-1).name), iflargefile);
frame_now  = read_B_snapshot(fullfile(path, B(iStart).name),   iflargefile);

%% ---------------- Main loop: centers k = iStart..iEnd ----------------
for k = iStart : iEnd

    fprintf('\n====================================================\n');
    fprintf('Center (%d/%d): %s\n', k, N_files, B(k).name);
    fprintf('Prev: %s | Next: %s\n', B(k-1).name, B(k+1).name);
    fprintf('====================================================\n');

    tic;

    frame_next = read_B_snapshot(fullfile(path, B(k+1).name), iflargefile);

    % Combine 3 frames
    [Xp, Yp, Zp, Up, Vp, Wp, axp, ayp, azp, T, weight, ...
        UNCUp, UNCVp, UNCWp, UNCAxp, UNCAyp, UNCAzp] = ...
        combine_three_frames(frame_last, frame_now, frame_next);

    % SVD interpolation
    [svd_v, svd_a] = run_svd_on_grid( ...
        Xg, Yg, Zg, Region, scale, ...
        Xp, Yp, Zp, Up, Vp, Wp, axp, ayp, azp, T, ...
        delta, nu, scale_unit, Nm_v, Nm_a, weight, ...
        UNCUp, UNCVp, UNCWp, cutoff_smooth);

    % CCM step
    [ccm_v, ccm_a] = run_ccm_step( ...
        Xg, Yg, Zg, Region, scale, ...
        Xp, Yp, Zp, Up, Vp, Wp, axp, ayp, azp, T, ...
        UNCUp, UNCVp, UNCWp, UNCAxp, UNCAyp, UNCAzp, ...
        svd_v, svd_a, delta, nu, scale_unit, Nm_v, Nm_a, weight);

    ccm_a(isnan(ccm_a)) = 0;

    % Output index: use filename number if available, else file index
    stepIndex = Bnum(k);
    if isnan(stepIndex), stepIndex = k; end

    [~,btag,~] = fileparts(B(k).name);
    tag = [btag '_'];

    % Write CCM
%     write_vel_file( ...
%         fullfile(dir_save_ccm, [fileName_prefix_CCM tag 'CCM_Vel_', num2str(stepIndex,'%05d'), '.dat']), ...
%         plotline1, plotline2v, plotline3, ccm_v(:,1:16));

    write_vel_file( ...
        fullfile(dir_save_ccm, ['CCM_Vel', num2str(stepIndex,'%05d'), '.dat']), ...
        plotline1, plotline2v, plotline3, ccm_v(:,1:16));


%     write_acc_file( ...
%         fullfile(dir_save_ccm, [fileName_prefix_CCM tag 'CCM_Acc_', num2str(stepIndex,'%05d'), '.dat']), ...
%         plotline1, plotline2a, plotline3, ccm_a(:,1:6), totalN);



    % Write SVD (optional)
    if ifsave_svd
%         write_vel_file( ...
%             fullfile(dir_save_svd, [fileName_prefix_SVD tag 'Vel_', num2str(stepIndex,'%05d'), '.dat']), ...
%             plotline1, plotline2v, plotline3, svd_v(:,1:16));
% 
%         write_acc_file( ...
%             fullfile(dir_save_svd, [fileName_prefix_SVD tag 'Acc_', num2str(stepIndex,'%05d'), '.dat']), ...
%             plotline1, plotline2a, plotline3, svd_a(:,1:6), totalN);

        write_vel_file( ...
            fullfile(dir_save_svd, ['SVD_Vel_', num2str(stepIndex,'%05d'), '.dat']), ...
            plotline1, plotline2v, plotline3, svd_v(:,1:16));

%         write_acc_file( ...
%             fullfile(dir_save_svd, ['Acc_', num2str(stepIndex,'%05d'), '.dat']), ...
%             plotline1, plotline2a, plotline3, svd_a(:,1:6), totalN);

    end

    fprintf('Done %s in %.2f s\n', B(k).name, toc);

    % Slide window
    frame_last = frame_now;
    frame_now  = frame_next;
end

fprintf('\nAll done.\n');



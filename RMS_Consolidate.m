% =========================================================================
% SCRIPT: RMS_Consolidate.m
% DESCRIPTION: Joins the per-topology RMS files written by
%              Motif_RMS_vs_ks_AllTopologies.m,
%                  mat/RMS_vs_ks_<Osc>_<topology><time_tag>.mat   (13 files),
%              into one consolidated file,
%                  mat/RMS_SyncIntervals_<Osc>_AllTopologies<time_tag>.mat
%              (RMS_data, RMS_array, SyncResults, SyncSummary,
%              intervals_all, ...), and prints the general table
%              (topologies x schemes) with the same cleaning as
%              RMS_PlotAndTable.m: only finite intervals and, with
%              drop_point_intervals = true, no single-point intervals [a, a],
%              and no interval narrower than min_width (intervals_all keeps
%              the raw intervals). No figures or LaTeX tables:
%              RMS_PlotAndTable.m makes them from the consolidated file.
%              - time_tag: '' for 'original', '_T<t_end>' otherwise.
%              - All files must share the same ks grid (kf, N_ks); a file
%                with another grid is skipped with a warning.
%              - Missing topologies are listed. With allow_incomplete =
%                true their rows stay NaN / "--" (topo_done = false); with
%                false nothing is written until all 13 exist.
%              Run it on its own (set OsID and time_profile below), or it
%              is called at the end of Motif_RMS_vs_ks_AllTopologies.m,
%              which passes consolidate_OsID and consolidate_time_profile.
% =========================================================================
if exist('consolidate_OsID', 'var')
    OsID         = consolidate_OsID;             % Called from Motif_RMS_vs_ks_AllTopologies.m
    time_profile = consolidate_time_profile;
else
    clear; clc; close all;
    % =====================================================================
    % CONFIGURATION (stand-alone use)
    % =====================================================================
    OsID         = 1;                   % 1=Rossler, 2=Lorenz
    time_profile = 'medium';            % 'original' | 'short' | 'medium' | 'long'
end
allow_incomplete = false;               % true: write even if some topologies are missing
mat_folder       = 'mat';
% General table: same cleaning as RMS_PlotAndTable.m (finite ends; with
% drop_point_intervals = true, single-point intervals [a, a] are removed;
% intervals narrower than min_width are removed)
drop_point_intervals = true;
min_width            = 0.1;             % minimum width of a valid interval (ks units; 0 = no limit)
num_fmt              = '%.2f';          % Number format of interval ends

OscType = {'Rossler', 'Lorenz', 'Chen', 'HR', 'Duffing', 'vanderPol'};
FMtype  = OscType{OsID};                % Same names as OscillatorType.m
switch lower(time_profile)
    case 'original', time_tag = '';
    case 'short',    time_tag = '_T1000';
    case 'medium',   time_tag = '_T2500';
    case 'long',     time_tag = '_T5000';
    otherwise, error('time_profile must be ''original'', ''short'', ''medium'' or ''long''.');
end

topology_list = { ...
    'chain', 'auxiliar', 'relay', 'mixed_master', 'competitive', 'mixed_slave', ...
    'ring_bidir_chain', 'ring_bidir_double_mixed', 'ring_bidir_mixed', ...
    'ring_unidir_chain', 'ring_bidir_auxiliar', 'ring_unidir_auxiliar', ...
    'ring_bidir_competitive'};
motif_names = {'Linear chain', 'Linear auxiliary', 'Linear relay', ...
    'Linear mixed master', 'Linear competitive', 'Linear mixed slave', ...
    'Ring bidirectional chain', 'Ring bidirectional double mixed', ...
    'Ring bidirectional mixed', 'Ring unidirectional chain', ...
    'Ring bidirectional auxiliary', 'Ring unidirectional auxiliary', ...
    'Ring bidirectional competitive'};
num_topo   = numel(topology_list);

% =========================================================================
% READ THE PER-TOPOLOGY FILES
% =========================================================================
RMS_data      = struct();
RMS_array     = [];
ks_vec        = [];
SyncResults   = struct('topology', {}, 'motif_name', {}, 'edges', {}, 'rms_e', {}, ...
                       'sync_intervals', {}, 'SyncTable', {});
intervals_all = cell(num_topo, 9);
interval_txt  = repmat("--", num_topo, 9);
topo_done     = false(1, num_topo);
run_time_min  = NaN(1, num_topo);
date_computed = strings(1, num_topo);
for topology_idx = 1:num_topo
    fname = fullfile(mat_folder, sprintf('RMS_vs_ks_%s_%s%s.mat', FMtype, topology_list{topology_idx}, time_tag));
    if ~isfile(fname)
        continue;
    end
    T = load(fname);
    if isempty(ks_vec)                  % The first file found sets the grid
        ks_vec       = T.ks_vec;
        kf           = T.kf;
        N_ks         = T.N_ks;
        sync_tol     = T.sync_tol;
        tspan        = T.tspan;
        t_transient  = T.t_transient;
        scheme_names = T.scheme_names;
        OscPar       = T.OscPar;
        RMS_array    = NaN(num_topo, 9, N_ks);
    elseif T.N_ks ~= N_ks || T.kf ~= kf
        warning('%s uses kf = %g, N_ks = %d (expected kf = %g, N_ks = %d): skipped.', ...
            fname, T.kf, T.N_ks, kf, N_ks);
        continue;
    end
    RMS_array(topology_idx, :, :) = reshape(T.RMS_topo, 1, 9, N_ks);
    RMS_data.(topology_list{topology_idx}) = T.RMS_data;
    SyncResults(topology_idx).topology       = T.topology_name;
    SyncResults(topology_idx).motif_name     = T.motif_name;
    SyncResults(topology_idx).edges          = T.edges;
    SyncResults(topology_idx).rms_e          = T.RMS_topo;
    SyncResults(topology_idx).sync_intervals = T.sync_intervals;
    SyncResults(topology_idx).SyncTable      = T.SyncTable;
    intervals_all(topology_idx, :) = T.sync_intervals(:).';
    interval_txt(topology_idx, :)  = T.interval_col(:).';
    run_time_min(topology_idx)     = T.run_time_min;
    date_computed(topology_idx)    = string(T.date_computed);
    topo_done(topology_idx)        = true;
    fprintf('Read: %s (%s, %.1f min)\n', fname, T.date_computed, T.run_time_min);
end

% =========================================================================
% CHECK
% =========================================================================
if ~any(topo_done)
    error('No file mat/RMS_vs_ks_%s_<topology>%s.mat found.', FMtype, time_tag);
end
if ~all(topo_done)
    fprintf('Missing topologies: %s\n', strjoin(topology_list(~topo_done), ', '));
    if ~allow_incomplete
        warning('%d of %d topologies missing: consolidated file NOT written (set allow_incomplete = true to write it anyway).', ...
            nnz(~topo_done), num_topo);
        return;
    end
end

% =========================================================================
% GENERAL TABLE (Command Window) AND CONSOLIDATED .MAT
% =========================================================================
% Valid intervals (as in RMS_PlotAndTable.m); intervals_all keeps the raw ones
valid_intervals = cell(num_topo, 9);
for topology_idx = find(topo_done)
    for CasoBs = 1:9
        I = intervals_all{topology_idx, CasoBs};
        if ~isempty(I)
            I = I(all(isfinite(I), 2), :);                 % Drop NaN/Inf rows
            if drop_point_intervals
                I = I(I(:, 2) > I(:, 1), :);               % Drop [a, a]
            end
            I = I(I(:, 2) - I(:, 1) >= min_width - 1e-9, :); % Drop narrow pieces
            interval_txt(topology_idx, CasoBs) = intervalText(I, num_fmt);
        end
        valid_intervals{topology_idx, CasoBs} = I;
    end
end
has_valid = ~cellfun(@isempty, valid_intervals);
SyncSummary = array2table(interval_txt, 'VariableNames', cellstr(scheme_names), ...
                          'RowNames', motif_names);
fprintf('\nGeneral table - synchronization intervals of ks (RMS < %g), %s, time profile ''%s'':\n', ...
    sync_tol, FMtype, time_profile);
disp(SyncSummary(topo_done, :));
fprintf('Schemes with synchronization   : %s\n', strjoin(cellstr(scheme_names(any(has_valid, 1))), ', '));
fprintf('Schemes without synchronization: %s\n\n', strjoin(cellstr(scheme_names(~any(has_valid, 1))), ', '));

filename_all = fullfile(mat_folder, sprintf('RMS_SyncIntervals_%s_AllTopologies%s.mat', FMtype, time_tag));
save(filename_all, 'RMS_data', 'RMS_array', 'scheme_names', ...
    'SyncResults', 'SyncSummary', 'intervals_all', ...
    'ks_vec', 'sync_tol', 'kf', 'N_ks', 'tspan', 't_transient', 'time_profile', ...
    'topology_list', 'motif_names', 'FMtype', 'OscPar', ...
    'topo_done', 'run_time_min', 'date_computed');
fprintf('Saved: %s (%d of %d topologies, total %.1f h)\n', filename_all, ...
    nnz(topo_done), num_topo, sum(run_time_min, 'omitnan') / 60);

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function str = intervalText(I, fmt)
    % Console text: "[a, b]"; several separated by "; " (as RMS_PlotAndTable.m)
    if isempty(I)
        str = "--";
        return;
    end
    parts = strings(size(I, 1), 1);
    for r = 1:size(I, 1)
        parts(r) = sprintf(['[' fmt ', ' fmt ']'], I(r, 1), I(r, 2));
    end
    str = join(parts, "; ");
end

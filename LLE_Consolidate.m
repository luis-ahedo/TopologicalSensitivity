% =========================================================================
% SCRIPT: LLE_Consolidate.m
% DESCRIPTION: Joins the per-topology LLE files written by
%              MotifLLE_vs_ks_AllTopologies.m,
%                  mat/LLE_vs_ks_<Osc>_<topology>.mat   (13 files),
%              into one consolidated file,
%                  mat/LLE_vs_ks_<Osc>_AllTopologies.mat,
%              which LLE_PlotAndTable.m reads to make the figures and the
%              LaTeX table.
%              - All files must share the same ks grid (kf, N_ks); a file
%                with another grid is skipped with a warning.
%              - Prints the table of synchronization intervals of ks
%                (LLE < lle_tol) of the topologies read, with the same
%                rules as LLE_PlotAndTable.m.
%              - Missing topologies are listed. With allow_incomplete =
%                true their rows stay NaN (topo_done = false); with false
%                nothing is saved until all 13 exist.
%              Run it on its own (set OsID below), or it is called at the
%              end of MotifLLE_vs_ks_AllTopologies.m (consolidate_at_end),
%              which passes its oscillator in consolidate_OsID.
% =========================================================================
if exist('consolidate_OsID', 'var')
    OsID = consolidate_OsID;            % Called from MotifLLE_vs_ks_AllTopologies.m
else
    clear; clc; close all;
    % =====================================================================
    % CONFIGURATION (stand-alone use)
    % =====================================================================
    OsID = 2;                           % 1=Rossler, 2=Lorenz
end 
allow_incomplete = false;               % true: save even if some topologies are missing
mat_folder       = 'mat';
% Console table of the synchronization intervals (same rules as LLE_PlotAndTable.m)
lle_tol              = 0;               % Synchronized where LLE < lle_tol
interp_crossings     = true;            % Linear interpolation of the zero crossings
drop_point_intervals = true;            % Drop intervals made of a single ks point
num_fmt              = '%.2f';          % Number format of interval ends

OscType = {'Rossler', 'Lorenz', 'Chen', 'HR', 'Duffing', 'vanderPol'};
FMtype  = OscType{OsID};                % Same names as OscillatorType.m

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
num_topo = numel(topology_list);

% =========================================================================
% READ THE PER-TOPOLOGY FILES
% =========================================================================
LLE_data      = struct();
LLE_array     = [];
ks_vec        = [];
topo_done     = false(1, num_topo);
run_time_min  = NaN(1, num_topo);
date_computed = strings(1, num_topo);
edges_all     = cell(1, num_topo);
LLE_coarse_array   = [];            % Coarse sweep (before edge refinement)
refined_mask_array = [];            % true where the ks point was recomputed with opts_refine
opts_refine_all    = cell(1, num_topo);
for topology_idx = 1:num_topo
    fname = fullfile(mat_folder, sprintf('LLE_vs_ks_%s_%s.mat', FMtype, topology_list{topology_idx}));
    if ~isfile(fname)
        continue;
    end
    T = load(fname);
    if isempty(ks_vec)                  % The first file found sets the grid
        ks_vec       = T.ks_vec;
        kf           = T.kf;
        N_ks         = T.N_ks;
        scheme_names = T.scheme_names;
        osc_names    = T.osc_names;
        OscPar       = T.OscPar;
        opts         = T.opts;
        LLE_array    = NaN(num_topo, 9, N_ks);
        LLE_coarse_array   = NaN(num_topo, 9, N_ks);
        refined_mask_array = false(num_topo, 9, N_ks);
    elseif T.N_ks ~= N_ks || T.kf ~= kf
        warning('%s uses kf = %g, N_ks = %d (expected kf = %g, N_ks = %d): skipped.', ...
            fname, T.kf, T.N_ks, kf, N_ks);
        continue;
    end
    LLE_array(topology_idx, :, :) = reshape(T.LLE_topo, 1, 9, N_ks);
    LLE_data.(topology_list{topology_idx}) = T.LLE_data;
    edges_all{topology_idx}     = T.edges;
    run_time_min(topology_idx)  = T.run_time_min;
    date_computed(topology_idx) = string(T.date_computed);
    topo_done(topology_idx)     = true;
    % Edge-refinement record (files from older versions may lack it)
    if isfield(T, 'LLE_coarse')
        LLE_coarse_array(topology_idx, :, :) = reshape(T.LLE_coarse, 1, 9, N_ks);
    end
    if isfield(T, 'refined_mask')
        refined_mask_array(topology_idx, :, :) = reshape(T.refined_mask, 1, 9, N_ks);
    end
    if isfield(T, 'opts_refine')
        opts_refine_all{topology_idx} = T.opts_refine;
    end
    fprintf('Read: %s (%s, %.1f min, %d refined ks points)\n', fname, T.date_computed, ...
        T.run_time_min, nnz(refined_mask_array(topology_idx, :, :)));
end

% =========================================================================
% CHECK AND SAVE
% =========================================================================
if ~any(topo_done)
    error('No file mat/LLE_vs_ks_%s_<topology>.mat found.', FMtype);
end

% Table of the synchronization intervals of the topologies read
% (CasoBs = 3*(i-1) + j for b_i -> c_j)
scheme_txt = strings(1, 9);
for CasoBs = 1:9
    scheme_txt(CasoBs) = sprintf('b%d->c%d', ceil(CasoBs / 3), mod(CasoBs - 1, 3) + 1);
end
SyncSummary = array2table(repmat("--", num_topo, 9), ...
    'VariableNames', cellstr(strrep(scheme_txt, '->', '_to_')), 'RowNames', motif_names);
has_valid = false(num_topo, 9);
for topology_idx = find(topo_done)
    for CasoBs = 1:9
        I = syncIntervalsLLE(ks_vec, squeeze(LLE_array(topology_idx, CasoBs, :)).', lle_tol, ...
            interp_crossings, drop_point_intervals);
        has_valid(topology_idx, CasoBs) = ~isempty(I);
        SyncSummary{topology_idx, CasoBs} = intervalText(I, num_fmt);
    end
end
fprintf('\nSynchronization intervals of ks (LLE < %g), %s:\n', lle_tol, FMtype);
disp(SyncSummary(topo_done, :));
fprintf('Schemes with an interval   : %s\n', strjoin(cellstr(scheme_txt(any(has_valid, 1))), ', '));
fprintf('Schemes without an interval: %s\n\n', strjoin(cellstr(scheme_txt(~any(has_valid, 1))), ', '));
if ~all(topo_done)
    fprintf('Missing topologies: %s\n', strjoin(topology_list(~topo_done), ', '));
    if ~allow_incomplete
        warning('%d of %d topologies missing: consolidated file NOT written (set allow_incomplete = true to write it anyway).', ...
            nnz(~topo_done), num_topo);
        return;
    end
end

filename_all = fullfile(mat_folder, sprintf('LLE_vs_ks_%s_AllTopologies.mat', FMtype));
save(filename_all, 'LLE_data', 'LLE_array', 'ks_vec', 'kf', 'N_ks', 'scheme_names', ...
    'topology_list', 'motif_names', 'osc_names', 'edges_all', 'FMtype', 'OsID', 'OscPar', ...
    'opts', 'run_time_min', 'date_computed', 'topo_done', ...
    'LLE_coarse_array', 'refined_mask_array', 'opts_refine_all');
fprintf('Saved: %s (%d of %d topologies, kf = %g, N_ks = %d, total %.1f h)\n', filename_all, ...
    nnz(topo_done), num_topo, kf, N_ks, sum(run_time_min, 'omitnan') / 60);

% =========================================================================
% LOCAL FUNCTIONS (same as LLE_PlotAndTable.m)
% =========================================================================
function I = syncIntervalsLLE(ks, L, tol, do_interp, drop_points)
    % Runs of ks where L < tol (NaN = not synchronized). Ends are
    % interpolated to the crossing L = tol when the neighbour is finite.
    % Right end = Inf when the run reaches the last ks point.
    mask  = L < tol;
    d     = diff([0, mask, 0]);
    i_on  = find(d == 1);
    i_off = find(d == -1) - 1;
    I = zeros(0, 2);
    N = numel(ks);
    for r = 1:numel(i_on)
        i1 = i_on(r);
        i2 = i_off(r);
        if drop_points && i1 == i2 && i2 < N
            continue;
        end
        a = ks(i1);
        if do_interp && i1 > 1 && isfinite(L(i1-1))
            a = ks(i1-1) + (tol - L(i1-1)) * (ks(i1) - ks(i1-1)) / (L(i1) - L(i1-1));
        end
        if i2 == N
            b = Inf;                                % Still synchronized at kf
        else
            b = ks(i2);
            if do_interp && isfinite(L(i2+1))
                b = ks(i2) + (tol - L(i2)) * (ks(i2+1) - ks(i2)) / (L(i2+1) - L(i2));
            end
        end
        I(end+1, :) = [a, b]; %#ok<AGROW>
    end
end

function str = intervalText(I, fmt)
    % Console text: "[a, b]" or "[a, Inf)"; several separated by "; "
    if isempty(I)
        str = "--";
        return;
    end
    parts = strings(size(I, 1), 1);
    for r = 1:size(I, 1)
        if isinf(I(r, 2))
            parts(r) = sprintf(['[' fmt ', Inf)'], I(r, 1));
        else
            parts(r) = sprintf(['[' fmt ', ' fmt ']'], I(r, 1), I(r, 2));
        end
    end
    str = join(parts, "; ");
end

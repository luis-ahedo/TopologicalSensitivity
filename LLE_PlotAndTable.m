% =========================================================================
% SCRIPT: LLE_PlotAndTable.m
% DESCRIPTION: Reads mat/LLE_vs_ks_<Osc>_AllTopologies.mat, made by
%              LLE_Consolidate.m from the per-topology files of
%              MotifLLE_vs_ks_AllTopologies.m (topologies not computed are
%              skipped and listed), for one oscillator and
%              1) plots LLE(ks) of the 9 coupling schemes in a 3x3
%                 tiledlayout for every topology, saving each figure to
%                 eps/ (figures are closed after saving);
%              2) detects the synchronization intervals (LLE < 0) and
%                 writes a LaTeX table: rows = topologies, columns = only
%                 the schemes with at least one interval.
%              Interval types: closed [a, b] when LLE returns above zero
%              inside [0, kf]; unbounded [a, \infty) when LLE is still
%              negative at kf (last point of the sweep).
% =========================================================================
clear; clc; close all;

% =========================================================================
% STEP 1: CONFIGURATION
% =========================================================================
OsID = 2;                           % 1=Rossler, 2=Lorenz
[~, ~, FMtype] = OscillatorType(OsID);   % Name used in the .mat file

mat_folder = 'mat';
eps_folder = 'eps';
tex_folder = 'tex';
filename_tex = fullfile(tex_folder, sprintf('LLE_SyncIntervals_%s_AllTopologies.tex', FMtype));

lle_tol              = 0;       % Synchronized where LLE < lle_tol
interp_crossings     = true;    % Linear interpolation of the zero crossings
drop_point_intervals = true;    % Drop intervals made of a single ks point
use_sideways         = true;    % true: sidewaystable, false: table [htbp]
num_fmt              = '%.2f';  % Number format of interval ends
col_sync    = [0.10 0.65 0.10]; % Curve colour where LLE < lle_tol (synchronization, green)
col_nosync  = [0.85 0.10 0.10]; % Curve colour where LLE >= lle_tol (no synchronization, red)

font_title     = 30;   % Global title (tiledlayout)
font_subtitle  = 30;   % Title of each subplot
font_xlabel    = 40;
font_ylabel    = 40;
font_ticks     = 20;

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
osc_names = {'R\"{o}ssler', 'Lorenz'};
num_topo  = numel(topology_list);
first_ring = 7;                     % Topologies 7..13 are rings

% =========================================================================
% STEP 2: LOAD THE CONSOLIDATED FILE (made by LLE_Consolidate.m)
%   LLE_array(topology, CasoBs, ks); rows of topologies not computed are
%   NaN and topo_done(topology) = false.
% =========================================================================
filename_mat = fullfile(mat_folder, sprintf('LLE_vs_ks_%s_AllTopologies.mat', FMtype));
if ~isfile(filename_mat)
    error('%s not found: run LLE_Consolidate.m first.', filename_mat);
end
S = load(filename_mat);
ks_vec    = S.ks_vec;
kf        = S.kf;
N_ks      = S.N_ks;
LLE_array = S.LLE_array;
if isfield(S, 'topo_done')
    topo_done = S.topo_done;
else
    topo_done = ~all(isnan(reshape(LLE_array, num_topo, [])), 2).';
end
fprintf('Loaded: %s (%d of %d topologies)\n', filename_mat, nnz(topo_done), num_topo);
if ~all(topo_done)
    fprintf('Not computed (skipped): %s\n', strjoin(topology_list(~topo_done), ', '));
end

scheme_tex = strings(1, 9);
scheme_txt = strings(1, 9);
for CasoBs = 1:9
    [bs, cs] = SISObc(CasoBs, 3);
    scheme_tex(CasoBs) = sprintf('$\\mathbf{b}_%d \\to \\mathbf{c}_%d$', find(bs), find(cs));
    scheme_txt(CasoBs) = sprintf('b%d->c%d', find(bs), find(cs));
end
osc_name = osc_names{OsID};

% =========================================================================
% STEP 3: ONE 3x3 FIGURE PER TOPOLOGY -> eps/ (figure closed after saving)
% =========================================================================
if ~exist(eps_folder, 'dir')
    mkdir(eps_folder);
end
for topology_idx = 1:num_topo
    if ~topo_done(topology_idx)
        continue;                               % Not computed yet
    end
    fig = figure('Visible', 'off');
    tl  = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for CasoBs = 1:9
        ax = nexttile(tl);
        % Green where LLE < lle_tol, red elsewhere; each colour also takes the
        % neighbouring point so the two pieces join at the crossings
        L  = reshape(LLE_array(topology_idx, CasoBs, :), 1, []);
        ms = L < lle_tol;
        mn = ~ms & ~isnan(L);
        Ls = L;  Ls(~(ms | [ms(2:end) false] | [false ms(1:end-1)])) = NaN;
        Ln = L;  Ln(~(mn | [mn(2:end) false] | [false mn(1:end-1)])) = NaN;
        plot(ax, ks_vec, Ln, '-', 'Color', col_nosync, 'LineWidth', 1.5);
        hold(ax, 'on');
        plot(ax, ks_vec, Ls, '-', 'Color', col_sync, 'LineWidth', 1.5);
        yline(ax, 0, 'k--', 'LineWidth', 1.0);
        grid(ax, 'on');
        xlim(ax, [0 kf]);
        set(ax, 'FontSize', font_ticks);
        title(ax, scheme_tex(CasoBs), 'Interpreter', 'latex', 'FontSize', font_subtitle);
    end
    xlabel(tl, '$k_s$', 'Interpreter', 'latex', 'FontSize', font_xlabel);
    ylabel(tl, '$\lambda^\perp_{\max}$', 'Interpreter', 'latex', 'FontSize', font_ylabel);
    title(tl, sprintf('%s - %s', osc_name, motif_names{topology_idx}), ...
        'Interpreter', 'latex', 'FontSize', font_title);

    filename_eps = fullfile(eps_folder, sprintf('Motif_LLE_vs_ks_AllSchemes_%s_%s_N%d.eps', ...
        FMtype, topology_list{topology_idx}, N_ks));
    exportgraphics(fig, filename_eps, 'ContentType', 'vector');
    close(fig);
    fprintf('Saved: %s\n', filename_eps);
end

% =========================================================================
% STEP 4: SYNCHRONIZATION INTERVALS (LLE < lle_tol)
%   intervals_all{topology, CasoBs} = N x 2 [a, b]; b = Inf if unbounded
% =========================================================================
intervals_all = cell(num_topo, 9);
for topology_idx = find(topo_done)
    for CasoBs = 1:9
        L = squeeze(LLE_array(topology_idx, CasoBs, :)).';
        intervals_all{topology_idx, CasoBs} = syncIntervalsLLE(ks_vec, L, lle_tol, ...
            interp_crossings, drop_point_intervals);
    end
end

has_valid = ~cellfun(@isempty, intervals_all);
keep_cols = find(any(has_valid, 1));
drop_cols = setdiff(1:9, keep_cols);

% Console table (all schemes)
SyncSummary = array2table(repmat("--", num_topo, 9), ...
    'VariableNames', cellstr(strrep(scheme_txt, '->', '_to_')), ...
    'RowNames', motif_names);
for topology_idx = find(topo_done)
    for CasoBs = 1:9
        SyncSummary{topology_idx, CasoBs} = intervalText(intervals_all{topology_idx, CasoBs}, num_fmt);
    end
end
SyncSummary = SyncSummary(topo_done, :);        % Computed topologies only
fprintf('\nSynchronization intervals of ks (LLE < %g), %s:\n', lle_tol, FMtype);
disp(SyncSummary);
fprintf('Kept columns   : %s\n', strjoin(cellstr(scheme_txt(keep_cols)), ', '));
fprintf('Removed columns: %s\n', strjoin(cellstr(scheme_txt(drop_cols)), ', '));

% =========================================================================
% STEP 5: LATEX TABLE (only schemes with intervals)
% =========================================================================
if isempty(keep_cols)
    warning('No scheme has a synchronization interval. No table written.');
    return;
end
if ~exist(tex_folder, 'dir')
    mkdir(tex_folder);
end
if use_sideways
    env = 'sidewaystable';  env_opt = '';
else
    env = 'table';          env_opt = '[htbp]';
end

fid = fopen(filename_tex, 'w');
fprintf(fid, '%% Generated by LLE_PlotAndTable.m from %s (%d of %d topologies)\n', ...
    filename_mat, nnz(topo_done), num_topo);
fprintf(fid, '%% Requires booktabs; sidewaystable is provided by sn-jnl (or \\usepackage{rotating})\n');
fprintf(fid, '\\begin{%s}%s\n', env, env_opt);
fprintf(fid, '\\centering\n');
fprintf(fid, ['\\caption{Intervals of $k_s$ with $\\lambda_{\\max} < %g$ (stable complete synchronization) ', ...
    'for the %s oscillator, per motif and coupling scheme; $[a, b]$: bounded interval, ', ...
    '$[a, \\infty)$: $\\lambda_{\\max} < %g$ up to $k_s = %g$. Schemes without synchronization omitted.}\n'], ...
    lle_tol, osc_name, lle_tol, kf);
fprintf(fid, '\\label{tab:lle-sync-%s}\n', FMtype);
fprintf(fid, '\\small\n');
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('c', 1, numel(keep_cols)));
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Motif & %s \\\\\n', strjoin(cellstr(scheme_tex(keep_cols)), ' & '));
fprintf(fid, '\\midrule\n');
rows_tex   = find(topo_done);                 % Computed topologies only
first_ring_row = find(rows_tex >= first_ring, 1);
for topology_idx = rows_tex
    if any(rows_tex < first_ring) && ~isempty(first_ring_row) && topology_idx == rows_tex(first_ring_row)
        fprintf(fid, '\\midrule\n');            % Linear | Ring separator
    end
    cells = strings(1, numel(keep_cols));
    for c = 1:numel(keep_cols)
        cells(c) = intervalTex(intervals_all{topology_idx, keep_cols(c)}, num_fmt);
    end
    fprintf(fid, '%s & %s \\\\\n', motif_names{topology_idx}, strjoin(cellstr(cells), ' & '));
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{%s}\n', env);
fclose(fid);
fprintf('Saved: %s\n', filename_tex);

% =========================================================================
% LOCAL FUNCTIONS
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

function str = intervalTex(I, fmt)
    % LaTeX cell: "--", "$[a, b]$", "$[a, \infty)$", or several stacked
    if isempty(I)
        str = "--";
        return;
    end
    parts = strings(size(I, 1), 1);
    for r = 1:size(I, 1)
        if isinf(I(r, 2))
            parts(r) = sprintf(['$[' fmt ', \\infty)$'], I(r, 1));
        else
            parts(r) = sprintf(['$[' fmt ', ' fmt ']$'], I(r, 1), I(r, 2));
        end
    end
    if numel(parts) == 1
        str = parts;
    else
        str = "\begin{tabular}[c]{@{}c@{}}" + join(parts, " \\ ") + "\end{tabular}";
    end
end
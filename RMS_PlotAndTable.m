% =========================================================================
% SCRIPT: RMS_PlotAndTable.m
% DESCRIPTION: Reads mat/RMS_SyncIntervals_<Osc>_AllTopologies<time_tag>.mat,
%              made by RMS_Consolidate.m from the per-topology files of
%              Motif_RMS_vs_ks_AllTopologies.m (topologies not computed are
%              skipped and listed), for one oscillator and
%              1) plots rms(e)(ks) of the 9 coupling schemes in a 3x3
%                 tiledlayout for every topology, on a log scale with the
%                 sync_tol line and the diverged runs shaded, saving each
%                 figure to eps/ (figures are closed after saving);
%              2) prints the synchronization intervals (rms(e) < sync_tol)
%                 and writes the LaTeX table of the manuscript: rows =
%                 topologies, columns = only the schemes with at least one
%                 valid interval,
%                     tex/RMS_SyncIntervals_<Osc>_AllTopologies<time_tag>.tex
%              A valid interval has finite ends, positive width (with
%              drop_point_intervals = true) and is at least min_width wide.
% =========================================================================
clear; clc; close all;

% =========================================================================
% STEP 1: CONFIGURATION
% =========================================================================
OsID         = 1;                   % 1=Rossler, 2=Lorenz
time_profile = 'medium';            % Same profile as Motif_RMS_vs_ks_AllTopologies.m: 'original' | 'short' | 'medium' | 'long'

mat_folder = 'mat';
eps_folder = 'eps';
tex_folder = 'tex';

make_figures         = true;        % One 3x3 figure per topology
write_reduced_table  = true;        % tex/..._AllTopologies<time_tag>.tex (schemes with synchronization only)
drop_point_intervals = true;        % true: [a, a] (single ks point) is not a valid interval
min_width            = 0.1;         % minimum width of a valid interval (ks units; 0 = no limit)
use_sideways         = true;        % true: sidewaystable, false: table [htbp]
num_fmt              = '%.2f';      % Number format of interval ends
rms_floor            = 1e-16;       % rms(e) values below this (e.g. exactly 0, converged early) are drawn at it
rms_tick_step        = 4;           % Decades between y ticks of the log axis
line_width           = 3.5;         % Line width of the rms(e) curves
col_curve            = [0 0 0];            % Curve colour (black)
col_tol              = [0.85 0.10 0.10];   % Colour of the dashed sync_tol line (red)
col_div              = [1 0.85 0.85];      % Shading where the run diverged (rms(e) = NaN)

font_title     = 30;   % Global title (tiledlayout)
font_subtitle  = 30;   % Title of each subplot
font_xlabel    = 40;
font_ylabel    = 40;
font_ticks     = 20;

OscType = {'Rossler', 'Lorenz', 'Chen', 'HR', 'Duffing', 'vanderPol'};   % Same names as OscillatorType.m
FMtype  = OscType{OsID};
osc_tex = {'R\"{o}ssler', 'Lorenz', 'Chen', 'Hindmarsh-Rose', 'Duffing', 'van der Pol'};
switch lower(time_profile)
    case 'original', time_tag = '';
    case 'short',    time_tag = '_T1000';
    case 'medium',   time_tag = '_T2500';
    case 'long',     time_tag = '_T5000';
    otherwise, error('time_profile must be ''original'', ''short'', ''medium'' or ''long''.');
end
first_ring = 7;                     % Topologies 7..13 are rings

% =========================================================================
% STEP 2: LOAD THE CONSOLIDATED FILE (made by RMS_Consolidate.m)
%   RMS_array(topology, CasoBs, ks); intervals_all{topology, CasoBs};
%   rows of topologies not computed are NaN / empty, topo_done = false.
% =========================================================================
filename_mat = fullfile(mat_folder, sprintf('RMS_SyncIntervals_%s_AllTopologies%s.mat', FMtype, time_tag));
if ~isfile(filename_mat)
    error('%s not found: run RMS_Consolidate.m first.', filename_mat);
end
S = load(filename_mat);
ks_vec        = S.ks_vec;
kf            = S.kf;
sync_tol      = S.sync_tol;
RMS_array     = S.RMS_array;
intervals_all = S.intervals_all;
topology_list = S.topology_list;
motif_names   = S.motif_names;
scheme_names  = S.scheme_names;
topo_done     = S.topo_done;
num_topo      = numel(topology_list);
fprintf('Loaded: %s (%d of %d topologies)\n', filename_mat, nnz(topo_done), num_topo);
if ~all(topo_done)
    fprintf('Not computed (skipped): %s\n', strjoin(topology_list(~topo_done), ', '));
end

scheme_tex = strings(1, 9);
for CasoBs = 1:9
    bc = sscanf(char(scheme_names(CasoBs)), 'b%d_to_c%d');          % [i; j] of b_i -> c_j
    scheme_tex(CasoBs) = sprintf('$\\mathbf{b}_%d \\to \\mathbf{c}_%d$', bc(1), bc(2));
end
rows_tex = find(topo_done);                                          % computed topologies
first_ring_row = rows_tex(find(rows_tex >= first_ring, 1));          % empty if no ring was computed

% =========================================================================
% STEP 3: ONE 3x3 FIGURE PER TOPOLOGY -> eps/ (figure closed after saving)
% =========================================================================
if make_figures
    if ~exist(eps_folder, 'dir'), mkdir(eps_folder); end
    dks = ks_vec(2) - ks_vec(1);
    for topology_idx = rows_tex
        % Same log y-range in every tile of this topology
        Rt   = RMS_array(topology_idx, :, :);
        rfin = Rt(isfinite(Rt));
        yr   = [rms_floor, 10^ceil(log10(max([rfin(:); 1])))];
        ye   = log10(yr(1)):rms_tick_step:log10(yr(2));
        fig = figure('Visible', 'off');
        tl  = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        for CasoBs = 1:9
            ax = nexttile(tl);
            hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
            set(ax, 'YScale', 'log');
            Rv  = reshape(RMS_array(topology_idx, CasoBs, :), 1, []);
            div = isnan(Rv);                                % Diverged (blow-up) runs
            if any(div)                                     % Shaded: diverged, no rms(e)
                d  = diff([0, div, 0]);
                on = find(d == 1);  off = find(d == -1) - 1;
                for q = 1:numel(on)
                    x1 = ks_vec(on(q)) - dks / 2;  x2 = ks_vec(off(q)) + dks / 2;
                    patch(ax, [x1 x2 x2 x1], yr([1 1 2 2]), col_div, 'EdgeColor', 'none');
                end
            end
            Rp = max(Rv, rms_floor);  Rp(div) = NaN;        % max(NaN, floor) would hide divergence
            plot(ax, ks_vec, Rp, '-', 'Color', col_curve, 'LineWidth', line_width);
            yline(ax, sync_tol, '--', 'Color', col_tol, 'LineWidth', 1.5);
            xlim(ax, [0 kf]);  ylim(ax, yr);
            set(ax, 'YTick', 10.^ye, 'YTickLabel', arrayfun(@(e) sprintf('$10^{%d}$', e), ye, ...
                'UniformOutput', false), 'Layer', 'top', 'TickLabelInterpreter', 'latex', ...
                'FontSize', font_ticks);
            title(ax, scheme_tex(CasoBs), 'Interpreter', 'latex', 'FontSize', font_subtitle);
        end
        xlabel(tl, '$k_s$', 'Interpreter', 'latex', 'FontSize', font_xlabel);
        ylabel(tl, '$RMS(\mathbf{e})$', 'Interpreter', 'latex', 'FontSize', font_ylabel);
        title(tl, sprintf('%s - %s', osc_tex{OsID}, motif_names{topology_idx}), ...
            'Interpreter', 'latex', 'FontSize', font_title);
        filename_eps = fullfile(eps_folder, sprintf('Motif_RMS_vs_ks_%s_%s%s.eps', ...
            FMtype, topology_list{topology_idx}, time_tag));
        exportgraphics(fig, filename_eps, 'ContentType', 'vector');
        close(fig);
        fprintf('Saved: %s\n', filename_eps);
    end
end

% =========================================================================
% STEP 4: SYNCHRONIZATION INTERVALS (rms(e) < sync_tol)
%   valid_intervals{topology, CasoBs}: finite ends, positive width
% =========================================================================
valid_intervals = cell(num_topo, 9);
for topology_idx = rows_tex
    for CasoBs = 1:9
        I = intervals_all{topology_idx, CasoBs};
        if ~isempty(I)
            I = I(all(isfinite(I), 2), :);                 % Drop NaN/Inf rows
            if drop_point_intervals
                I = I(I(:, 2) > I(:, 1), :);               % Drop [a, a]
            end
            I = I(I(:, 2) - I(:, 1) >= min_width - 1e-9, :); % Drop narrow pieces
        end
        valid_intervals{topology_idx, CasoBs} = I;
    end
end
has_valid = ~cellfun(@isempty, valid_intervals);
keep_cols = find(any(has_valid, 1));
drop_cols = setdiff(1:9, keep_cols);

% Console table (all schemes)
SyncSummary = array2table(repmat("--", num_topo, 9), 'VariableNames', cellstr(scheme_names), ...
    'RowNames', motif_names);
for topology_idx = rows_tex
    for CasoBs = 1:9
        SyncSummary{topology_idx, CasoBs} = intervalText(valid_intervals{topology_idx, CasoBs}, num_fmt);
    end
end
fprintf('\nSynchronization intervals of ks (rms(e) < %g), %s, time profile ''%s'':\n', sync_tol, FMtype, time_profile);
disp(SyncSummary(topo_done, :));
fprintf('Schemes with synchronization   : %s\n', strjoin(scheme_names(keep_cols), ', '));
fprintf('Schemes without synchronization: %s\n', strjoin(scheme_names(drop_cols), ', '));

% =========================================================================
% STEP 5: LATEX TABLE (schemes with synchronization only)
% =========================================================================
if ~exist(tex_folder, 'dir'), mkdir(tex_folder); end
if write_reduced_table
    if isempty(keep_cols)
        warning('No scheme has a valid synchronization interval: reduced table not written.');
    else
        filename_red = fullfile(tex_folder, sprintf('RMS_SyncIntervals_%s_AllTopologies%s.tex', FMtype, time_tag));
        writeTable(filename_red, filename_mat, keep_cols, valid_intervals, rows_tex, first_ring, first_ring_row, ...
            scheme_tex, motif_names, num_fmt, use_sideways, ...
            sprintf(['Intervals of $k_s$ with $rms(\\mathbf{e}) < %g$ for the %s oscillator, per motif and ', ...
            'coupling scheme (schemes without synchronization omitted).'], sync_tol, osc_tex{OsID}), ...
            sprintf('tab:sync-%s-reduced', FMtype), num_topo);
    end
end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function writeTable(filename_tex, filename_mat, cols, intervals, rows_tex, first_ring, first_ring_row, ...
                    scheme_tex, motif_names, num_fmt, use_sideways, caption, label, num_topo)
    % LaTeX table: rows = computed topologies, columns = schemes in cols
    if use_sideways
        env = 'sidewaystable';  env_opt = '';
    else
        env = 'table';          env_opt = '[htbp]';
    end
    fid = fopen(filename_tex, 'w');
    fprintf(fid, '%% Generated by RMS_PlotAndTable.m from %s (%d of %d topologies)\n', ...
        filename_mat, numel(rows_tex), num_topo);
    fprintf(fid, '%% Requires booktabs; sidewaystable is provided by sn-jnl (or \\usepackage{rotating})\n');
    fprintf(fid, '\\begin{%s}%s\n', env, env_opt);
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{%s}\n', caption);
    fprintf(fid, '\\label{%s}\n', label);
    fprintf(fid, '\\small\n');
    fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('c', 1, numel(cols)));
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'Motif & %s \\\\\n', strjoin(cellstr(scheme_tex(cols)), ' & '));
    fprintf(fid, '\\midrule\n');
    for topology_idx = rows_tex
        if ~isempty(first_ring_row) && topology_idx == first_ring_row && any(rows_tex < first_ring)
            fprintf(fid, '\\midrule\n');            % Linear | Ring separator
        end
        cells = strings(1, numel(cols));
        for c = 1:numel(cols)
            cells(c) = intervalTex(intervals{topology_idx, cols(c)}, num_fmt);
        end
        fprintf(fid, '%s & %s \\\\\n', motif_names{topology_idx}, strjoin(cellstr(cells), ' & '));
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{%s}\n', env);
    fclose(fid);
    fprintf('Saved: %s\n', filename_tex);
end

function str = intervalText(I, fmt)
    % Console text: "[a, b]"; several separated by "; "
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

function str = intervalTex(I, fmt)
    % LaTeX cell: "--", "$[a, b]$", or intervals stacked vertically
    if isempty(I)
        str = "--";
        return;
    end
    parts = strings(size(I, 1), 1);
    for r = 1:size(I, 1)
        parts(r) = sprintf(['$[' fmt ', ' fmt ']$'], I(r, 1), I(r, 2));
    end
    if numel(parts) == 1
        str = parts;
    else
        str = "\begin{tabular}[c]{@{}c@{}}" + join(parts, " \\ ") + "\end{tabular}";
    end
end
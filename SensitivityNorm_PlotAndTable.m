% =========================================================================
% SCRIPT: SensitivityNorm_PlotAndTable.m
% DESCRIPTION: Reads mat/SensitivityNorm_vs_ks_<Osc>_AllTopologies_N<N>_T<T_end>.mat,
%              made by SensitivityNorm_Consolidate.m from the per-case files
%              of SensitivityNorm_vs_ks_AllTopologies.m, for one oscillator:
%              1) prints, per coupling scheme with at least one computed
%                 case, the topological sensitivity index <S_q> (dB) of every
%                 topology with its uncertainty and swept ks interval;
%              2) draws one bar chart per scheme: one group of six horizontal
%                 bars (da2 db2 dc2 da3 db3 dc3) per topology, groups
%                 separated by thin lines, saved to eps/ (figures are closed
%                 after saving; all share the same x-axis, common_xlim);
%              3) writes ONE LaTeX table with all schemes: rows = topology x
%                 scheme, columns = <S_q> of the six mismatches (optionally
%                 the swept ks interval, show_interval):
%                 tex/TopologicalSensitivity_<Osc>_AllSchemes_N<N>_T<T_end>.tex
%              <S_q> = mean over the valid ks (internal RMS check passed and
%              convergence test met) of 20 log10 of the step-1 curve
%              sqrt((1/3) sum_ij S_q(e_ij)^2), as computed by
%              SensitivityNorm_Consolidate.m (TS_tab, TS_unc_tab, TS_gap_tab).
%              Topologies without a result for a scheme (no valid interval,
%              all ks flagged) are left out of that figure and table and
%              listed in the command window.
% =========================================================================
clear; clc; close all;

% =========================================================================
% STEP 1: CONFIGURATION
% =========================================================================
OsID         = 1;               % 1=Rossler, 2=Lorenz
N_sweep      = 25;              % Same N_sweep as the computation
time_profile = 'medium';        % Same time_profile as the computation: 'short' | 'medium' | 'long'
scheme_run   = [];              % [] = every scheme with a result; or CasoBs values, e.g. [1 5]

mat_folder = 'mat';
eps_folder = 'eps';
tex_folder = 'tex';

make_figures   = true;          % One bar chart per scheme
write_tables   = true;          % One LaTeX table with all schemes
show_interval  = false;         % true: column with the swept ks interval of each case
common_xlim    = true;          % true: same x-axis in every figure
show_title     = true;          % false for the paper (the caption describes the scheme)
show_values    = false;         % value of <S_q> at the end of each bar (crowded with many topologies)
show_errorbars = false;         % +/- 2 x uncertainty (usually smaller than the line width)
show_unc_table = false;         % true: cells "value +/- uncertainty"; false: values, largest uncertainty in the caption
use_sideways   = false;         % true: sidewaystable, false: table [htbp]
num_fmt        = '%.1f';        % Number format of <S_q> (dB)
int_fmt        = '%.2f';        % Number format of the swept interval ends

fig_width   = 560;              % Figure size (pixels): width, height per topology group and fixed part
fig_h_group = 48;
fig_h_fixed = 150;
font_title  = 16;
font_xlabel = 16;
font_legend = 13;
font_ticks  = 13;
font_values = 9;

OscType   = {'Rossler', 'Lorenz', 'Chen', 'HR', 'Duffing', 'vanderPol'};   % Same names as OscillatorType.m
FMtype    = OscType{OsID};
osc_tex   = {'R\"{o}ssler', 'Lorenz', 'Chen', 'Hindmarsh-Rose', 'Duffing', 'van der Pol'};
topo_short = {'Lin. chain', 'Lin. aux.', 'Lin. relay', 'Lin. mix. mast.', 'Lin. comp.', ...
    'Lin. mix. slav.', 'Ring bidir. chain', 'Ring bidir. dbl. mix.', 'Ring bidir. mix.', ...
    'Ring unidir. chain', 'Ring bidir. aux.', 'Ring unidir. aux.', 'Ring bidir. comp.'};
par_tex   = {'\delta a_2', '\delta b_2', '\delta c_2', '\delta a_3', '\delta b_3', '\delta c_3'};
par_txt   = {'da2', 'db2', 'dc2', 'da3', 'db3', 'dc3'};
par_tex_tab = {'\delta \mathcal{a}_2', '\delta \mathcal{b}_2', '\delta \mathcal{c}_2', ...
    '\delta \mathcal{a}_3', '\delta \mathcal{b}_3', '\delta \mathcal{c}_3'};   % LaTeX table (paper notation; MATLAB figures keep par_tex)
switch lower(time_profile)
    case 'short',  T_end = 1000;
    case 'medium', T_end = 2500;
    case 'long',   T_end = 5000;
    otherwise, error('time_profile must be ''short'', ''medium'' or ''long''.');
end
first_ring = 7;                 % Topologies 7..13 are rings

% =========================================================================
% STEP 2: LOAD THE UNIFIED FILE (made by SensitivityNorm_Consolidate.m)
%   TS_tab(topology, CasoBs, q), TS_unc_tab, TS_gap_tab, status_tab,
%   ks_a_tab, ks_b_tab (13 x 9 [x 6])
% =========================================================================
filename_mat = fullfile(mat_folder, sprintf('SensitivityNorm_vs_ks_%s_AllTopologies_N%d_T%d.mat', FMtype, N_sweep, T_end));
if ~isfile(filename_mat)
    error('%s not found: run SensitivityNorm_Consolidate.m first.', filename_mat);
end
S = load(filename_mat, 'TS_tab', 'TS_unc_tab', 'TS_gap_tab', 'status_tab', 'ks_a_tab', 'ks_b_tab', ...
    'topology_list', 'motif_names', 'scheme_bc', 'lle_file');
TS         = S.TS_tab;
TS_unc     = S.TS_unc_tab;
TS_gap     = S.TS_gap_tab;
status_tab = S.status_tab;
num_topo   = numel(S.topology_list);
fprintf('Loaded: %s\n', filename_mat);

has_result = any(isfinite(TS), 3);                          % 13 x 9
if isempty(scheme_run)
    scheme_run = find(any(has_result, 1));
end
if isempty(scheme_run)
    error('No case with a result in %s.', filename_mat);
end

% Intervals that reach kf (end of the LLE grid) are marked in the tables
kf = NaN;
if isfield(S, 'lle_file') && isfile(S.lle_file)
    Lk = load(S.lle_file, 'kf');
    kf = Lk.kf;
end

% =========================================================================
% STEP 3: COMMAND WINDOW TABLE
% =========================================================================
for CasoBs = scheme_run
    bc = S.scheme_bc(CasoBs, :);
    fprintf('\n<S_q> (dB) +/- uncertainty, %s, b%d -> c%d (CasoBs = %d)\n', FMtype, bc(1), bc(2), CasoBs);
    fprintf('%-22s %-16s', 'Topology', 'Swept ks');
    fprintf(' %13s', par_txt{:});
    fprintf('\n');
    for ti = 1:num_topo
        fprintf('%-22s', topo_short{ti});
        if ~has_result(ti, CasoBs)
            fprintf(' %s\n', status_tab(ti, CasoBs));
            continue;
        end
        fprintf(' %-16s', sprintf('[%.2f, %.2f]', S.ks_a_tab(ti, CasoBs), S.ks_b_tab(ti, CasoBs)));
        fprintf(' %6.2f +/-%4.2f', [reshape(TS(ti, CasoBs, :), 1, []); reshape(TS_unc(ti, CasoBs, :), 1, [])]);
        fprintf('\n');
    end
end

% =========================================================================
% STEP 4: BAR CHARTS, one per scheme -> eps/ (figure closed after saving)
% =========================================================================
base_col = [0.000 0.447 0.741;                      % a
            0.850 0.325 0.098;                      % b
            0.494 0.184 0.556];                     % c
col_q = [base_col; 0.45 * base_col + 0.55];         % node 2 (dark), node 3 (light)
leg_q = cellfun(@(s) ['$' s '$'], par_tex, 'UniformOutput', false);
v_all = reshape(TS(:, scheme_run, :), [], 1);
v_all = v_all(isfinite(v_all));
x_common = barLimits(v_all, show_values);          % margins: no bar or label on the frame

if make_figures
    if ~exist(eps_folder, 'dir'), mkdir(eps_folder); end
    for CasoBs = scheme_run
        bc    = S.scheme_bc(CasoBs, :);
        shown = find(has_result(:, CasoBs)).';
        if isempty(shown), continue; end
        nV = numel(shown);
        Y  = reshape(TS(shown, CasoBs, :), nV, 6);
        if common_xlim
            x_lim = x_common;
        else
            x_lim = barLimits(Y(:), show_values);
        end

        fig = figure('Visible', 'off', 'Units', 'pixels', ...
            'Position', [100, 100, fig_width, fig_h_fixed + fig_h_group * nV]);
        ax = axes(fig);
        hold(ax, 'on');
        hb = barh(ax, 1:nV+1, [Y; NaN(1, 6)], 'grouped', 'BarWidth', 0.85, ...   % NaN row: grouped bars also for one topology
            'EdgeColor', 'k', 'LineWidth', 0.4);
        for q = 1:6
            hb(q).FaceColor = col_q(q, :);
        end
        for q = 1:6
            pos = hb(q).XEndPoints(1:nV);
            val = Y(:, q).';
            if show_errorbars
                errorbar(ax, val, pos, 2 * reshape(TS_unc(shown, CasoBs, q), 1, []), 'horizontal', ...
                    'Color', 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 2, 'HandleVisibility', 'off');
            end
            if show_values
                for iv = 1:nV
                    if ~isfinite(val(iv)), continue; end
                    if val(iv) >= 0                 % label right of the bar end
                        x_txt = val(iv) + 0.4;  h_al = 'left';
                    else                            % negative bar: label left of its end
                        x_txt = val(iv) - 0.4;  h_al = 'right';
                    end
                    text(ax, x_txt, pos(iv), sprintf(num_fmt, val(iv)), 'FontSize', font_values, ...
                        'Interpreter', 'latex', 'HorizontalAlignment', h_al, 'VerticalAlignment', 'middle');
                end
            end
        end
        for iv = 1:nV-1                                 % separators between topology groups
            yline(ax, iv + 0.5, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
        xline(ax, 0, 'k-', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        grid(ax, 'on');
        ax.YGrid = 'off';
        set(ax, 'YDir', 'reverse', 'YTick', 1:nV, 'YTickLabel', topo_short(shown), ...
            'FontSize', font_ticks, 'TickLabelInterpreter', 'latex', 'Layer', 'top', 'Box', 'on');
        ylim(ax, [0.5, nV + 0.5]);
        xlim(ax, x_lim);
        xlabel(ax, '$\langle \bar{S}_q(\mathbf{e}) \rangle$ (dB)', 'Interpreter', 'latex', 'FontSize', font_xlabel);
        if show_title
            title(ax, sprintf('Sensitivity by topology, $\\mathbf{b}_%d \\to \\mathbf{c}_%d$', bc(1), bc(2)), ...
                'Interpreter', 'latex', 'FontSize', font_title);
        end
        legend(ax, hb, leg_q, 'Interpreter', 'latex', 'FontSize', font_legend, ...
            'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 6);
        ax.Toolbar.Visible = 'off';
        filename_eps = fullfile(eps_folder, sprintf('TopologicalSensitivity_AllTopologies_%s_b%dc%d.eps', ...
            FMtype, bc(1), bc(2)));
        exportgraphics(fig, filename_eps, 'ContentType', 'vector');
        close(fig);
        fprintf('Saved: %s\n', filename_eps);
    end
end

% =========================================================================
% STEP 5: ONE LATEX TABLE, all schemes
%   rows = topology x scheme (topologies in order, schemes of scheme_run
%   with a result), columns = [swept ks interval] + <S_q> of the six
%   mismatches
% =========================================================================
if write_tables
    if ~exist(tex_folder, 'dir'), mkdir(tex_folder); end
    if use_sideways
        env = 'sidewaystable';  env_opt = '';
    else
        env = 'table';          env_opt = '[htbp]';
    end
    rows_t = find(any(has_result(:, scheme_run), 2)).';
    first_ring_row = rows_t(find(rows_t >= first_ring, 1));
    unc_max = max(reshape(TS_unc(rows_t, scheme_run, :), [], 1), [], 'omitnan');
    gap_max = max(reshape(TS_gap(rows_t, scheme_run, :), [], 1), [], 'omitnan');

    filename_tex = fullfile(tex_folder, sprintf('TopologicalSensitivity_%s_AllSchemes_N%d_T%d.tex', ...
        FMtype, N_sweep, T_end));
    fid = fopen(filename_tex, 'w');
    fprintf(fid, '%% Generated by SensitivityNorm_PlotAndTable.m from %s\n', filename_mat);
    fprintf(fid, '%% Requires booktabs\n');
    fprintf(fid, '\\begin{%s}%s\n', env, env_opt);
    fprintf(fid, '\\centering\n');
    if show_unc_table
        unc_txt = 'Each value $\\pm$ its Monte Carlo uncertainty.';
    else
        unc_txt = sprintf(['Monte Carlo uncertainty at most %.2f~dB; difference between the two ', ...
            'estimators of the index, the average of the per-$k_s$ medians ', ...
            '[Eq.~\\eqref{eq_Sq_global}] and the median over the realizations of the ', ...
            'per-realization index, at most %.2f~dB.'], unc_max, gap_max);
    end
    if show_interval
        int_txt = sprintf(' $k_s$ swept over the given interval ($^{*}$: interval reaching $k_f = %g$).', kf);
    else
        int_txt = '';
    end
    fprintf(fid, ['\\caption{Topological sensitivity index $\\langle \\bar{S}_q(\\mathbf{e}) \\rangle$ (dB) of the %s oscillator ', ...
        'per motif and coupling scheme.%s %s Abbreviations: Lin.\\ = linear; mix.\\ = mixed; mast.\\ = master; ', ...
        'slav.\\ = slave; bidir.\\ = bidirectional; unidir.\\ = unidirectional; dbl.\\ = double; ', ...
        'aux.\\ = auxiliary; comp.\\ = competitive.}\n'], osc_tex{OsID}, int_txt, unc_txt);
    fprintf(fid, '\\label{tab:topo-sens-%s}\n', FMtype);
    fprintf(fid, '\\footnotesize\n');
    par_head = strjoin(cellfun(@(s) ['$' s '$'], par_tex_tab, 'UniformOutput', false), ' & ');
    if show_interval
        fprintf(fid, '\\begin{tabular}{lllcccccc}\n\\toprule\n');
        fprintf(fid, 'Motif & Scheme & $k_s$ & %s \\\\\n', par_head);
    else
        fprintf(fid, '\\begin{tabular}{llcccccc}\n\\toprule\n');
        fprintf(fid, 'Motif & Scheme & %s \\\\\n', par_head);
    end
    fprintf(fid, '\\midrule\n');
    for ti = rows_t
        if ~isempty(first_ring_row) && ti == first_ring_row && any(rows_t < first_ring)
            fprintf(fid, '\\midrule\n');                % Linear | Ring separator
        end
        for CasoBs = scheme_run
            if ~has_result(ti, CasoBs), continue; end
            bc = S.scheme_bc(CasoBs, :);
            cells = strings(1, 6);
            for q = 1:6
                v = TS(ti, CasoBs, q);
                if ~isfinite(v)
                    cells(q) = "--";
                elseif show_unc_table
                    cells(q) = sprintf(['$' num_fmt ' \\pm ' num_fmt '$'], v, TS_unc(ti, CasoBs, q));
                else
                    cells(q) = sprintf(['$' num_fmt '$'], v);
                end
            end
            sch_txt = sprintf('$\\mathbf{b}_%d \\to \\mathbf{c}_%d$', bc(1), bc(2));
            if show_interval
                ks_txt = sprintf(['$[' int_fmt ', ' int_fmt ']$'], S.ks_a_tab(ti, CasoBs), S.ks_b_tab(ti, CasoBs));
                if isfinite(kf) && abs(S.ks_b_tab(ti, CasoBs) - kf) < 1e-9
                    ks_txt = [ks_txt '$^{*}$']; %#ok<AGROW>
                end
                fprintf(fid, '%s & %s & %s & %s \\\\\n', topo_short{ti}, sch_txt, ks_txt, strjoin(cellstr(cells), ' & '));
            else
                fprintf(fid, '%s & %s & %s \\\\\n', topo_short{ti}, sch_txt, strjoin(cellstr(cells), ' & '));
            end
        end
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{%s}\n', env);
    fclose(fid);
    fprintf('Saved: %s\n', filename_tex);
end

% -------------------------------------------------------------------------
% x-axis of a bar chart of <S_q>: the zero line is included, the ends are
% multiples of 5 dB and a margin is left on each side (wider when the
% values are written at the bar ends), so no bar or label touches the frame
% -------------------------------------------------------------------------
function x_lim = barLimits(v, show_values)
    v     = v(isfinite(v));
    v_lo  = min([v(:); 0]);
    v_hi  = max([v(:); 0]);
    pad   = (0.04 + 0.10 * show_values) * max(v_hi - v_lo, 1);
    x_lim = [floor((v_lo - pad * (v_lo < 0)) / 5) * 5, ceil((v_hi + pad * (v_hi > 0)) / 5) * 5];
    if x_lim(2) <= x_lim(1), x_lim(2) = x_lim(1) + 5; end
end

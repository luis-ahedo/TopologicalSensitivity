% =========================================================================
% SCRIPT: PlotValidationRMS.m
% DESCRIPTION: Redraws the validation figure (measured vs predicted RMS of
%              the synchronization errors, Fig. 7) from the files saved by
%              ValidateSensitivity_vs_ks.m, without integrating again:
%                  mat/ValidationRMS_vs_ks_<Osc>_<topology>_B<C>_<dq>_r<r>.mat
%              one file per mismatch size r (one run of
%              ValidateSensitivity_vs_ks.m per value of mismatch_ratio).
%              The measured RMS(e_ij) of every r is divided by r and drawn
%              against the single first-order prediction
%                  RMS(e_ij)/r = x_ref * 10^(S_q(e_ij)/20),
%              so the curves collapse onto it where the first-order theory
%              holds and fan out where it fails. With one value of r the
%              figure is the original one (RMS, not divided by r).
%              Same layout and style as STEP 5 of ValidateSensitivity_vs_ks.m.
%              Output: eps/ValidationRMS_vs_ks_<tag>.eps (overwritten)
% =========================================================================
clear; clc; close all;

% =========================================================================
% STEP 1: CONFIGURATION
% =========================================================================
file_base   = 'ValidationRMS_vs_ks_Rossler_chain_B1_da2';   % file name without _r<r>.mat
r_list      = [0.01 0.05 0.1 0.2];  % mismatch sizes to draw (missing files are skipped)

show_title  = true;             % false for the paper (the caption describes the case)
show_band   = false;            % +/- 2 x uncertainty band of the prediction
y_log       = true;             % Logarithmic y-axis
y_max       = [];               % Linear axis only: upper y-limit ([] = from the measured values)
fig_size    = [100 100 600 800];
font_title  = 18;
font_label  = 16;
font_legend = 11;
font_ticks  = 13;

osc_list    = {'Rossler', 'Lorenz'};
osc_tex     = {'R\"{o}ssler', 'Lorenz'};
topology_list = {'chain', 'auxiliar', 'relay', 'mixed_master', 'competitive', 'mixed_slave', ...
    'ring_bidir_chain', 'ring_bidir_double_mixed', 'ring_bidir_mixed', ...
    'ring_unidir_chain', 'ring_bidir_auxiliar', 'ring_unidir_auxiliar', 'ring_bidir_competitive'};
motif_names = {'Linear chain', 'Linear auxiliary', 'Linear relay', ...
    'Linear mixed master', 'Linear competitive', 'Linear mixed slave', ...
    'Ring bidirectional chain', 'Ring bidirectional double mixed', ...
    'Ring bidirectional mixed', 'Ring unidirectional chain', ...
    'Ring bidirectional auxiliary', 'Ring unidirectional auxiliary', ...
    'Ring bidirectional competitive'};
param_names = {'da2', 'db2', 'dc2', 'da3', 'db3', 'dc3'};
par_letter  = {'a', 'b', 'c'};

% =========================================================================
% STEP 2: LOAD (one file per r)
% =========================================================================
V_all = {};
r_ok  = [];
for r = r_list
    f = fullfile('mat', sprintf('%s_r%g.mat', file_base, r));
    if ~isfile(f)
        warning('%s not found: r = %g skipped (run ValidateSensitivity_vs_ks.m with mismatch_ratio = %g).', f, r, r);
        continue;
    end
    V_all{end + 1} = load(f);                       %#ok<SAGROW>
    r_ok(end + 1)  = abs(V_all{end}.mismatch_ratio); %#ok<SAGROW>
    fprintf('Loaded: %s\n', f);
end
if isempty(V_all), error('No validation file found for %s.', file_base); end
nR        = numel(V_all);
norm_r    = nR > 1;                                 % divide by r only when several r are drawn
V         = V_all{1};                               % case description and prediction from the first file
ks_prof   = V.ks_prof(:);
scale     = 1;
if norm_r, scale = r_ok(1); end
pred      = V.pred / scale;                         % RMS/r: the same for every r
pred_lo   = V.pred_lo / scale;
pred_hi   = V.pred_hi / scale;
no_infl   = V.no_infl;
OsID       = V.OsID;
CasoBs     = V.CasoBs;
FMtype     = osc_list{OsID};
topology   = topology_list{V.topology_idx};
i_b        = ceil(CasoBs / 3);                      % same convention as ValidateSensitivity_vs_ks.m
j_c        = mod(CasoBs - 1, 3) + 1;
q_idx      = mod(V.mismatch_idx - 1, 3) + 1;
node_mm    = 2 + (V.mismatch_idx > 3);
num_pairs  = 3;

% =========================================================================
% STEP 3: FIGURE - measured vs predicted, one panel per error vector
% =========================================================================
col_pred = [0.850 0.325 0.098];
col_meas = [0 0.447 0.741];
if norm_r                                           % sequential blues, light (small r) to dark (large r)
    col_r = interp1([0 1], [0.55 0.75 0.95; 0 0.25 0.55], linspace(0, 1, nR));
else
    col_r = col_meas;
end
mk_r     = {'o', 's', '^', 'd', 'v', 'p'};
err_tex  = {'12', '13', '23'};
r_lbl    = arrayfun(@(r) sprintf('Measured, $r = %g\\,\\%%$', 100 * r), r_ok, 'UniformOutput', false);
if ~norm_r, r_lbl = {'Measured'}; end
fig = figure('Name', 'Measured vs predicted RMS error', 'Position', fig_size);
tl  = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
h_leg = gobjects(1, nR + 1);                        % common legend (first panel with a prediction)
for p = 1:num_pairs
    ax = nexttile(tl);
    hold(ax, 'on');
    h = gobjects(1, nR + 1);
    y_all = [];
    if ~no_infl(p)
        ok = isfinite(pred(:, p));
        if show_band && any(ok)
            fill(ax, [ks_prof(ok); flipud(ks_prof(ok))], [pred_lo(ok, p); flipud(pred_hi(ok, p))], ...
                col_pred, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        h(nR + 1) = plot(ax, ks_prof, pred(:, p), '--', 'Color', col_pred, 'LineWidth', 2);
        y_all = pred(:, p);
    end
    for k = 1:nR
        ks_m = V_all{k}.ks_measure(:);
        y_m  = V_all{k}.meas(:, p);
        if norm_r, y_m = y_m / r_ok(k); end
        h(k) = plot(ax, ks_m, y_m, ['-' mk_r{mod(k - 1, numel(mk_r)) + 1}], 'Color', col_r(k, :), ...
            'MarkerFaceColor', col_r(k, :), 'MarkerSize', 4, 'LineWidth', 1.5);
        y_all = [y_all; y_m];                       %#ok<AGROW>
    end
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontSize', font_ticks, 'TickLabelInterpreter', 'latex');
    ks_max = max(cellfun(@(W) max(W.ks_measure), V_all));
    xlim(ax, [0, ceil(max([ks_max; ks_prof]))]);
    if y_log
        set(ax, 'YScale', 'log');                   % Whole decades covering every measured and predicted value
        y_all = y_all(isfinite(y_all) & y_all > 0);
        if ~isempty(y_all)
            ylim(ax, 10.^[floor(log10(min(y_all))), ceil(log10(max(y_all)))]);
        end
    elseif ~isempty(y_max)
        ylim(ax, [0, y_max]);
    else
        y_top = max(y_all(isfinite(y_all)));
        if ~isempty(y_top) && y_top > 0
            ylim(ax, [0, ceil(1.1 * y_top / 0.5) * 0.5]);
        end
    end
    if norm_r
        ylabel(ax, sprintf('$\\mathrm{RMS}(\\mathbf{e}_{%s})/r$', err_tex{p}), 'Interpreter', 'latex', 'FontSize', font_label);
    else
        ylabel(ax, sprintf('$\\mathrm{RMS}(\\mathbf{e}_{%s})$', err_tex{p}), 'Interpreter', 'latex', 'FontSize', font_label);
    end
    text(ax, 0.97, 0.90, sprintf('(%s)', char('a' + p - 1)), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'Interpreter', 'latex', 'FontSize', font_ticks);
    if no_infl(p)
        text(ax, 0.03, 0.90, 'no influence predicted', 'Units', 'normalized', ...
            'Interpreter', 'latex', 'FontSize', font_legend);
    elseif ~any(isgraphics(h_leg))
        h_leg = h;                                  % handles of the common legend
    end
end
xlabel(tl, '$k_s$', 'Interpreter', 'latex', 'FontSize', font_label);
if ~any(isgraphics(h_leg)), h_leg = [h(1:nR), gobjects(1)]; end   % every panel without prediction
ok_leg = isgraphics(h_leg);
leg_txt = [r_lbl, {'Predicted'}];
lgd = legend(h_leg(ok_leg), leg_txt(ok_leg), 'Interpreter', 'latex', 'FontSize', font_legend, ...
    'NumColumns', min(3, nnz(ok_leg)));             % one legend for the three panels, below them (rows of 3 fit the figure width)
lgd.Layout.Tile = 'south';
if show_title
    if norm_r
        r_txt = sprintf('$\\delta %s_%d = r\\,%s$', par_letter{q_idx}, node_mm, par_letter{q_idx});
    else
        r_txt = sprintf('$\\delta %s_%d = %g\\,%s$', par_letter{q_idx}, node_mm, r_ok(1), par_letter{q_idx});
    end
    title(tl, sprintf('%s - %s, $\\mathbf{b}_%d \\to \\mathbf{c}_%d$: %s', osc_tex{OsID}, ...
        motif_names{V.topology_idx}, i_b, j_c, r_txt), 'Interpreter', 'latex', 'FontSize', font_title);
end
all_ax = findall(fig, 'Type', 'axes');
for kk = 1:numel(all_ax)
    all_ax(kk).Toolbar.Visible = 'off';
end

% =========================================================================
% STEP 4: SAVE (eps/)
% =========================================================================
if ~exist('eps', 'dir'), mkdir('eps'); end
if norm_r                                           % e.g. ..._da2_r1_5_10_20pct
    file_tag = sprintf('%s_%s_B%d_%s_r%spct', FMtype, topology, CasoBs, param_names{V.mismatch_idx}, ...
        strjoin(arrayfun(@(r) strrep(sprintf('%g', 100 * r), '.', 'p'), r_ok, 'UniformOutput', false), '_'));
else
    file_tag = sprintf('%s_%s_B%d_%s_r%g', FMtype, topology, CasoBs, param_names{V.mismatch_idx}, r_ok(1));
end
exportgraphics(fig, fullfile('eps', sprintf('ValidationRMS_vs_ks_%s.eps', file_tag)), 'ContentType', 'vector');
fprintf('Saved: eps/ValidationRMS_vs_ks_%s.eps\n', file_tag);

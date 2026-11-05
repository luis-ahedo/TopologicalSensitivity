% =========================================================================
% SCRIPT: CompareMotifTopologiesErrors.m
% DESCRIPTION: Compares the synchronization errors of several motif
%              configurations (topology + coupling scheme + ks), e.g. the
%              2 x 2 design {topology A, topology B} x {scheme 1, scheme 2},
%              under the SAME single parameter mismatch, for one or several
%              mismatch levels, and compares the measured RMS errors with the
%              first-order prediction from the LOCAL sensitivity index at the
%              ks of each configuration:
%                  RMS(e_ij) ~ x_ref * 10^(S_q(e_ij)/20) * r,
%              with S_q(e_ij) interpolated (in dB) on the sweep of the
%              consolidated file. ks = NaN selects ks*, the ks of the sweep at
%              which lambda_max is minimum over the converged ks (Eq. eq_Sq_star,
%              same rule as SensitivityStar.m). Self-contained (vector fields and coupling written
%              here), same conventions as the current pipeline:
%                - node k: dx_k/dt = F(x_k; q + dq_k) + ks*sum_{s->k} B*C*(x_s - x_k),
%                  node 1 nominal, single mismatch dq = r*q on node 2 or 3;
%                - start ON the synchronization manifold: x1 = x2 = x3 =
%                  a point of the node-1 attractor (same point for all
%                  configurations), so the errors are produced only by the
%                  mismatch;
%                - RMS window [t_transient, T_end] sampled every dt_sample
%                  (same as the sensitivity profiles, 'medium': [500, 2500]);
%                  vector RMS(e_ij) = ( sum_xyz mean(e^2) )^(1/2).
%              One figure, tiled layout (rows = mismatch levels, columns =
%              e12, e13, e23): ||e_ij(t)|| of every configuration over a
%              plotted segment of the window; console table of measured and
%              predicted RMS(e_ij) and their ratio; LaTeX table
%              tex/CompareErrorsTable_<tag>.tex and mat/CompareErrors_<tag>.mat.
% =========================================================================
clear; clc; close all;

% =========================================================================
% STEP 1: CONFIGURATION
% =========================================================================
OsID = 1;                       % 1 = Rossler, 2 = Lorenz

%   1 chain, 2 auxiliar, 3 relay, 4 mixed_master, 5 competitive,
%   6 mixed_slave, 7 ring_bidir_chain, 8 ring_bidir_double_mixed,
%   9 ring_bidir_mixed, 10 ring_unidir_chain, 11 ring_bidir_auxiliar,
%   12 ring_unidir_auxiliar, 13 ring_bidir_competitive
% One row per configuration: [topology, CasoBs, ks]; CasoBs = 3*(i-1) + j for b_i -> c_j;
% ks = NaN -> ks* of that configuration (needs the consolidated file)
configs = [1 1 NaN;             % linear chain, b1 -> c1
           1 5 NaN;             % linear chain, b2 -> c2
           8 1 NaN;             % ring bidirectional double mixed, b1 -> c1
           8 5 NaN];            % ring bidirectional double mixed, b2 -> c2

mismatch_idx    = 3;                     % 1..6 = da2, db2, dc2, da3, db3, dc3
mismatch_ratios = [0.05 0.10 0.20];      % r = delta q / q, one row of the figure per value

T_end       = 2500;             % Same window as the sensitivity profiles ('medium')
t_transient = 500;
dt_sample   = 0.01;
plot_window = [t_transient, t_transient + 200];   % segment of ||e(t)|| drawn in the figure
plot_dt     = 0.05;             % time step of the PLOTTED curves (the RMS uses dt_sample); keeps the vector EPS light
plot_pairs  = [1 2];            % plotted errors: 1 = e12, 2 = e13 (e23 = e13 - e12 is redundant); the RMS table keeps all three
common_ylim = true;             % true: same y-axis in all rows (levels comparable); false: one y-axis per row
ode_RelTol  = 1e-8;
ode_AbsTol  = 1e-10;
pool_T_transient = 500;         % node-1 starting point: point of the attractor after this transient
base_seed   = 42;

N_sweep      = 25;              % to read ks*, S_q(e_ij) and <S_q> from mat/SensitivityNorm_vs_ks_<Osc>_AllTopologies_N<N>_T<T>.mat
font_title  = 18;
font_label  = 16;
font_legend = 12;
font_ticks  = 13;
fig_size    = [100 100 900 400];
write_table = true;             % LaTeX table of measured vs predicted RMS -> tex/

% =========================================================================
% STEP 2: OSCILLATOR, TOPOLOGIES AND COUPLING (no external functions)
% =========================================================================
osc_list = {'Rossler', 'Lorenz'};
osc_tex  = {'R\"{o}ssler', 'Lorenz'};
FMtype   = osc_list{OsID};
switch FMtype
    case 'Rossler'
        OscPar = [0.2, 0.2, 5.7];
        Fnode  = @(x, p) [-x(2) - x(3); x(1) + p(1)*x(2); p(2) + x(3)*(x(1) - p(3))];
    case 'Lorenz'
        OscPar = [10, 28, 2];   % same values as OscillatorType.m
        Fnode  = @(x, p) [p(1)*(x(2) - x(1)); x(1)*(p(2) - x(3)) - x(2); x(1)*x(2) - p(3)*x(3)];
end
n = 3;
topology_list = {'chain', 'auxiliar', 'relay', 'mixed_master', 'competitive', 'mixed_slave', ...
    'ring_bidir_chain', 'ring_bidir_double_mixed', 'ring_bidir_mixed', ...
    'ring_unidir_chain', 'ring_bidir_auxiliar', 'ring_unidir_auxiliar', 'ring_bidir_competitive'};
topo_short = {'Lin. chain', 'Lin. aux.', 'Lin. relay', 'Lin. mix. mast.', 'Lin. comp.', ...
    'Lin. mix. slav.', 'Ring bidir. chain', 'Ring bidir. dbl. mix.', 'Ring bidir. mix.', ...
    'Ring unidir. chain', 'Ring bidir. aux.', 'Ring unidir. aux.', 'Ring bidir. comp.'};
edges_all = {[1 2; 2 3], [1 2; 1 3], [1 2; 2 1; 2 3; 3 2], [1 2; 2 3; 3 2], [1 2; 3 2], ...
    [2 1; 2 3; 3 2], [1 2; 2 1; 2 3; 3 2; 3 1; 1 3], [1 2; 2 3; 3 2; 3 1; 1 3], ...
    [1 2; 2 3; 3 2; 3 1], [1 2; 2 3; 3 1], [1 2; 1 3; 2 3; 3 2], [1 2; 1 3; 2 3], ...
    [1 2; 2 1; 2 3; 1 3]};
param_names = {'da2', 'db2', 'dc2', 'da3', 'db3', 'dc3'};
pair_lbl    = {'12', '13', '23'};
par_letter  = {'a', 'b', 'c'};
q_idx   = mod(mismatch_idx - 1, 3) + 1;
node_mm = 2 + (mismatch_idx > 3);

nC  = size(configs, 1);
cfg = struct('topo', num2cell(configs(:, 1)).', 'CasoBs', num2cell(configs(:, 2)).', 'ks', num2cell(configs(:, 3)).');

% Consolidated sensitivity file: ks*, local indices and x_ref
cons_file = fullfile('mat', sprintf('SensitivityNorm_vs_ks_%s_AllTopologies_N%d_T%d.mat', FMtype, N_sweep, T_end));
has_cons  = isfile(cons_file);
xref      = NaN;
if has_cons
    C    = load(cons_file);
    if ~isfield(C, 'xref')                          % norm_dB must be the normalized index
        error('%s has no xref: rerun SensitivityNorm_Consolidate.m on the shifted files.', cons_file);
    end
    xref = C.xref;
else
    warning('%s not found: no ks*, no prediction and no global index.', cons_file);
end
for c = 1:nC
    if isnan(cfg(c).ks)
        if ~has_cons, error('ks = NaN needs %s.', cons_file); end
        Rk = C.SN(cfg(c).topo, cfg(c).CasoBs);
        lam = Rk.lle_sweep(:);                      % same rule as SensitivityStar.m:
        lam(~(isfinite(lam) & any(isfinite(Rk.TS_curve_dB), 2))) = NaN;   % only ks with a converged index
        [~, i0] = min(lam);
        cfg(c).ks = Rk.ks_sweep(i0);                % ks* (Table tab:sens-star)
    end
end
for c = 1:nC
    i_b = ceil(cfg(c).CasoBs / 3);
    j_c = mod(cfg(c).CasoBs - 1, 3) + 1;
    BC  = zeros(n);  BC(i_b, j_c) = 1;              % B*C, same convention as SISObc.m
    cfg(c).BC    = BC;
    cfg(c).bc    = [i_b, j_c];
    cfg(c).edges = edges_all{cfg(c).topo};
    cfg(c).label = sprintf('%s ($\\mathbf{b}_%d \\to \\mathbf{c}_%d$, $k_s = %.2f$)', topo_short{cfg(c).topo}, i_b, j_c, cfg(c).ks);
end

% Node-1 starting point on its attractor (same for all configurations)
rng(base_seed, 'twister');
[~, Xp] = ode45(@(t, x) Fnode(x, OscPar), [0, pool_T_transient], rand(n, 1), ...
    odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
x1_0 = Xp(end, :).';
X0   = [x1_0; x1_0; x1_0];                          % on the synchronization manifold

t_samp = t_transient + (0:floor((T_end - t_transient) / dt_sample)) * dt_sample;
ode_opts = odeset('RelTol', ode_RelTol, 'AbsTol', ode_AbsTol, 'Events', @(t, X) blowupEvent(X, 1e3));

% Local indices S_q(e_ij) at the ks of each configuration (prediction) and global <S_q>
S_loc = NaN(nC, 3);                                 % dB, per error pair
S_glb = NaN(nC, 1);                                 % <S_q> (dB)
if has_cons
    for c = 1:nC
        Rk = C.SN(cfg(c).topo, cfg(c).CasoBs);
        if isempty(Rk.ks_sweep) || cfg(c).ks < min(Rk.ks_sweep) || cfg(c).ks > max(Rk.ks_sweep)
            warning('Configuration %d: ks = %g outside its sweep, no prediction.', c, cfg(c).ks);
        else
            for p = 1:3
                v = (mismatch_idx - 1) * 3 + p;
                S_loc(c, p) = interp1(Rk.ks_sweep(:), Rk.norm_dB(:, v), cfg(c).ks);
            end
        end
        S_glb(c) = C.TS_tab(cfg(c).topo, cfg(c).CasoBs, mismatch_idx);
        fprintf('Config %d %s B%d: ks = %.3f, S_q(e12, e13, e23) = %.1f, %.1f, %.1f dB, <S_q> = %.1f dB\n', ...
            c, topology_list{cfg(c).topo}, cfg(c).CasoBs, cfg(c).ks, S_loc(c, :), S_glb(c));
    end
end

% =========================================================================
% STEP 3: SIMULATIONS, all mismatch levels
% =========================================================================
nL      = numel(mismatch_ratios);
RMS_all = NaN(nL, nC, 3);                           % level x configuration x pair
PRED    = NaN(nL, nC, 3);                           % predicted RMS (local index)
E_all   = cell(nL, nC, 3);                           % ||e_ij(t)|| in the RMS window
fprintf('\n%-8s %-4s %10s %10s %10s\n', 'r', 'cfg', 'RMS(e12)', 'RMS(e13)', 'RMS(e23)');
for ir = 1:nL
    r  = mismatch_ratios(ir);
    DP = zeros(3, 3);
    DP(node_mm, q_idx) = r * OscPar(q_idx);
    P_all = repmat(OscPar, 3, 1) + DP;
    for c = 1:nC
        PRED(ir, c, :) = xref * 10.^(S_loc(c, :) / 20) * r;
        f_net = @(t, X) networkRHS(X, Fnode, P_all, cfg(c).edges, cfg(c).BC, cfg(c).ks, n);
        [t_out, X_out, te] = ode45(f_net, [0, t_samp], X0, ode_opts);
        if ~isempty(te) || numel(t_out) < numel(t_samp) + 1
            warning('Configuration %d diverged for r = %g.', c, r);
            continue;
        end
        Xs = X_out(2:end, :);
        e  = {Xs(:, 1:3) - Xs(:, 4:6), Xs(:, 1:3) - Xs(:, 7:9), Xs(:, 4:6) - Xs(:, 7:9)};
        for p = 1:3
            E_all{ir, c, p}   = vecnorm(e{p}, 2, 2);
            RMS_all(ir, c, p) = sqrt(sum(mean(e{p}.^2, 1)));
        end
        fprintf('%-8g %-4d %10.4f %10.4f %10.4f   predicted %8.4f %8.4f %8.4f   meas/pred %5.2f %5.2f %5.2f\n', ...
            r, c, squeeze(RMS_all(ir, c, :)), squeeze(PRED(ir, c, :)), squeeze(RMS_all(ir, c, :) ./ PRED(ir, c, :)));
    end
end

% =========================================================================
% STEP 4: ONE FIGURE, tiled layout nL x 3
%   rows = mismatch levels, columns = e12, e13, e23; ||e_ij(t)|| of both
%   configurations over plot_window; same y-axis in all rows (common_ylim)
%   or within each row
% =========================================================================
if ~exist('eps', 'dir'), mkdir('eps'); end
t_plot = t_samp(2:end).';
in_win = t_plot >= plot_window(1) & t_plot <= plot_window(2);
stride = max(1, round(plot_dt / dt_sample));
k_win  = find(in_win);
in_win = false(size(in_win));
in_win(k_win(1:stride:end)) = true;              % plotted points only (every plot_dt)
fig = figure('Name', 'Error comparison', 'Position', [fig_size(1:3), 280 * nL + 80]);
nP  = numel(plot_pairs);
tl  = tiledlayout(fig, nL, nP, 'TileSpacing', 'compact', 'Padding', 'compact');
hC = gobjects(1, nC);
cols = lines(nC);
styles = {'-', '--', '-.', ':'};
% y-range from the maximum of the PLOTTED errors (e23 excluded when not plotted)
E_plot = E_all(:, :, plot_pairs);
y_all  = max([0, cellfun(@(v) max([v(in_win); 0]), E_plot(~cellfun(@isempty, E_plot))).']);
for ir = 1:nL
    Er = reshape(E_plot(ir, :, :), nC, nP);
    y_max = max([0, cellfun(@(v) max([v(in_win); 0]), Er(~cellfun(@isempty, Er))).']);
    if common_ylim, y_max = y_all; end
    if y_max == 0, y_max = 1; end
    for p = plot_pairs
        ax = nexttile(tl);
        hold(ax, 'on');
        for c = 1:nC
            if ~isempty(E_all{ir, c, p})
                hC(c) = plot(ax, t_plot(in_win), E_all{ir, c, p}(in_win), styles{mod(c - 1, 4) + 1}, ...
                    'LineWidth', 1.2, 'Color', cols(c, :));
            end
        end
        grid(ax, 'on');  box(ax, 'on');
        set(ax, 'FontSize', font_ticks, 'TickLabelInterpreter', 'latex');
        xlim(ax, plot_window);
        ylim(ax, [0, 1.05 * y_max]);                % 5 % margin above the largest plotted error
        if ir == 1
            title(ax, sprintf('$\\mathbf{e}_{%s}$', pair_lbl{p}), 'Interpreter', 'latex', 'FontSize', font_label);
        end
        if ir == nL
            xlabel(ax, '$t$', 'Interpreter', 'latex', 'FontSize', font_label);
        else
            xticklabels(ax, {});
        end
        if p == plot_pairs(1)                       % mismatch level of the row
            ylabel(ax, {sprintf('$\\delta %s_%d = %g\\,%s$', par_letter{q_idx}, node_mm, mismatch_ratios(ir), par_letter{q_idx}), ...
                sprintf('$\\|\\mathbf{e}_{%s}(t)\\|$', pair_lbl{p})}, 'Interpreter', 'latex', 'FontSize', font_label);
        else
            ylabel(ax, sprintf('$\\|\\mathbf{e}_{%s}(t)\\|$', pair_lbl{p}), 'Interpreter', 'latex', 'FontSize', font_label);
        end
    end
end
ok_h = isgraphics(hC);
lgd = legend(hC(ok_h), {cfg(ok_h).label}, 'Interpreter', 'latex', 'FontSize', font_legend, 'NumColumns', 2);
lgd.Layout.Tile = 'south';
title(tl, sprintf('%s: mismatch $\\delta %s_%d$', osc_tex{OsID}, par_letter{q_idx}, node_mm), ...
    'Interpreter', 'latex', 'FontSize', font_title);
all_ax = findall(fig, 'Type', 'axes');
for kk = 1:numel(all_ax), all_ax(kk).Toolbar.Visible = 'off'; end
lvl_txt = strjoin(arrayfun(@(v) sprintf('%g', 100 * v), mismatch_ratios, 'UniformOutput', false), '_');
cfg_txt = strjoin(arrayfun(@(c) sprintf('%s_B%d', topology_list{cfg(c).topo}, cfg(c).CasoBs), 1:nC, 'UniformOutput', false), '_');
out_eps = fullfile('eps', sprintf('CompareErrors_%s_%s_%s_%spct.eps', FMtype, cfg_txt, param_names{mismatch_idx}, lvl_txt));
warning('off', 'MATLAB:print:ContentTypeImageSuggested');
exportgraphics(fig, out_eps, 'ContentType', 'vector');
fprintf('Saved: %s\n', out_eps);

% =========================================================================
% STEP 5: LATEX TABLE, measured vs predicted RMS(e_ij) -> tex/  and mat/
% =========================================================================
if write_table && has_cons
    if ~exist('tex', 'dir'), mkdir('tex'); end
    out_tex = fullfile('tex', sprintf('CompareErrorsTable_%s_%s_%s_%spct.tex', FMtype, cfg_txt, param_names{mismatch_idx}, lvl_txt));
    fid = fopen(out_tex, 'w');
    fprintf(fid, '%% Generated by CompareMotifTopologiesErrors.m from %s\n', cons_file);
    fprintf(fid, '\\begin{tabular}{lcccccccc}\n\\toprule\n');
    fprintf(fid, ['Motif & Scheme & $k_s$ & $S_q(\\mathbf{e}_{12})$, $S_q(\\mathbf{e}_{13})$ (dB) & $r$ & ', ...
        '\\multicolumn{2}{c}{$\\operatorname{RMS}(\\mathbf{e}_{12})$} & \\multicolumn{2}{c}{$\\operatorname{RMS}(\\mathbf{e}_{13})$} \\\\\n']);
    fprintf(fid, ' & & & & & meas. & pred. & meas. & pred. \\\\\n\\midrule\n');
    for c = 1:nC
        for ir = 1:nL
            if ir == 1
                head = sprintf('%s & $\\mathbf{b}_%d \\to \\mathbf{c}_%d$ & $%.2f$ & $%.1f$, $%.1f$', topo_short{cfg(c).topo}, ...
                    cfg(c).bc(1), cfg(c).bc(2), cfg(c).ks, S_loc(c, 1), S_loc(c, 2));
            else
                head = ' & & &';
            end
            fprintf(fid, '%s & $%g$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\\\\n', head, mismatch_ratios(ir), ...
                RMS_all(ir, c, 1), PRED(ir, c, 1), RMS_all(ir, c, 2), PRED(ir, c, 2));
        end
        if c < nC, fprintf(fid, '\\midrule\n'); end
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
    fclose(fid);
    fprintf('Saved: %s\n', out_tex);
end
if ~exist('mat', 'dir'), mkdir('mat'); end
save(fullfile('mat', sprintf('CompareErrors_%s_%s_%s_%spct.mat', FMtype, cfg_txt, param_names{mismatch_idx}, lvl_txt)), ...
    'configs', 'cfg', 'mismatch_idx', 'mismatch_ratios', 'RMS_all', 'PRED', 'S_loc', 'S_glb', 'xref', 'x1_0', ...
    'T_end', 't_transient', 'dt_sample', 'ode_RelTol', 'ode_AbsTol');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function dX = networkRHS(X, Fnode, P_all, edges, BC, ks, n)
    dX = zeros(3*n, 1);
    for k = 1:3
        dX((k-1)*n + (1:n)) = Fnode(X((k-1)*n + (1:n)), P_all(k, :));
    end
    for e = 1:size(edges, 1)
        s = edges(e, 1);
        k = edges(e, 2);
        dX((k-1)*n + (1:n)) = dX((k-1)*n + (1:n)) + ks * BC * (X((s-1)*n + (1:n)) - X((k-1)*n + (1:n)));
    end
end

function [value, isterminal, direction] = blowupEvent(X, bound)
    if any(~isfinite(X)), value = -1; else, value = bound - max(abs(X)); end
    isterminal = 1;
    direction  = 0;
end

% =========================================================================
% SCRIPT: ValidateSensitivity_vs_ks.m
% DESCRIPTION: Validation of the first-order sensitivity prediction
%              against the measured RMS synchronization error, vs ks.
%              Self-contained (vector fields and coupling written here);
%              it only READS the sensitivity profile saved by
%              SensitivityNorm_vs_ks.m:
%                mat/SensitivityNorm_vs_ks_<Osc>_<topology>_B<CasoBs>_N<N>_T<T>.mat
%
%   Predicted (Eq. eq_measured_predicte_normalized, normalized sensitivities):
%     RMS(e_ij) ~ x_ref * ||s_qj(e_ij)|| * |delta q_j / q_j|
%     ||s_qj(e_ij)|| = 10^(norm_dB / 20)   (linear, from the profile file)
%     The sensitivities are already scaled by the nominal q_j
%     and divided by the attractor scale x_ref (sbar = s * q_j / x_ref,
%     field xref of the profile file), so the factor is x_ref times the
%     RATIO r = delta q_j / q_j, NOT the absolute mismatch r * q_j.
%   Measured:
%     Full nonlinear network (3 nodes) with the single mismatch
%     delta q_j = r * q_j on node 2 or 3 (node 1 nominal), same topology,
%     coupling scheme, ks values, RMS window [t_transient, T_end] and
%     sampling step as the sensitivity profile. Each run starts ON the
%     synchronization manifold (x1 = x2 = x3 = point r of the node-1
%     attractor pool of the profile file), so the measured error is
%     produced only by the mismatch.
%     RMS(e_ij) = ( sum over x,y,z of mean(e_ij,comp^2) )^(1/2)  per run,
%     median over the n_runs runs (same statistic as the sensitivity).
% =========================================================================
clear; clc; close all;

% =========================================================================
% STEP 1: CONFIGURATION
% =========================================================================
OsID         = 1;               % 1 = Rossler, 2 = Lorenz
topology_idx = 1;               % 1 = chain (see topology_list below)
CasoBs       = 1;               % b_i -> c_j: CasoBs = 3*(i-1) + j
profile_file = '';              % '' = newest mat/SensitivityNorm_vs_ks_<Osc>_<topology>_B<CasoBs>_N*_T*.mat

mismatch_idx   = 1;             % 1..6 = da2, db2, dc2, da3, db3, dc3
mismatch_ratio = 0.1;           % r = delta q_j / q_j (10 %)

n_runs      = 10;               % Runs (pool points) per ks for the measured RMS
ks_measure  = [];               % [] = the ks of the profile; or a vector, e.g. linspace(0.2, 4, 60)
ode_RelTol  = 1e-8;
ode_AbsTol  = 1e-10;
err_bound   = 1e3;              % Blow-up bound on |x| (run excluded)

show_title  = true;             % false for the paper (the caption describes the case)
show_band   = false;            % +/- 2 x uncertainty band of the prediction (narrower than the line: off for the paper)
y_max       = [];               % Upper y-limit of every panel ([] = from the measured values of each panel; ignored if y_log)
y_log       = true;             % Logarithmic y-axis: interior agreement and edge growth both visible
fig_size    = [100 100 600 800];
font_title  = 18;
font_label  = 16;
font_legend = 11;
font_ticks  = 13;

% =========================================================================
% STEP 2: OSCILLATOR, TOPOLOGY AND COUPLING (no external functions)
% =========================================================================
osc_list  = {'Rossler', 'Lorenz'};
osc_tex   = {'R\"{o}ssler', 'Lorenz'};
FMtype    = osc_list{OsID};
switch FMtype
    case 'Rossler'
        OscPar = [0.2, 0.2, 5.7];
        Fnode  = @(x, p) [-x(2) - x(3); x(1) + p(1)*x(2); p(2) + x(3)*(x(1) - p(3))];
    case 'Lorenz'
        OscPar = [10, 28, 2];     % same values as OscillatorType.m
        Fnode  = @(x, p) [p(1)*(x(2) - x(1)); x(1)*(p(2) - x(3)) - x(2); x(1)*x(2) - p(3)*x(3)];
end
n = 3;

topology_list = {'chain', 'auxiliar', 'relay', 'mixed_master', 'competitive', 'mixed_slave', ...
    'ring_bidir_chain', 'ring_bidir_double_mixed', 'ring_bidir_mixed', ...
    'ring_unidir_chain', 'ring_bidir_auxiliar', 'ring_unidir_auxiliar', 'ring_bidir_competitive'};
motif_names = {'Linear chain', 'Linear auxiliary', 'Linear relay', ...
    'Linear mixed master', 'Linear competitive', 'Linear mixed slave', ...
    'Ring bidirectional chain', 'Ring bidirectional double mixed', ...
    'Ring bidirectional mixed', 'Ring unidirectional chain', ...
    'Ring bidirectional auxiliary', 'Ring unidirectional auxiliary', ...
    'Ring bidirectional competitive'};
edges_all = {[1 2; 2 3], [1 2; 1 3], [1 2; 2 1; 2 3; 3 2], [1 2; 2 3; 3 2], [1 2; 3 2], ...
    [2 1; 2 3; 3 2], [1 2; 2 1; 2 3; 3 2; 3 1; 1 3], [1 2; 2 3; 3 2; 3 1; 1 3], ...
    [1 2; 2 3; 3 2; 3 1], [1 2; 2 3; 3 1], [1 2; 1 3; 2 3; 3 2], [1 2; 1 3; 2 3], ...
    [1 2; 2 1; 2 3; 1 3]};
topology = topology_list{topology_idx};
edges    = edges_all{topology_idx};

% Same convention as SISObc.m: CasoBs -> B = e_i (input), C = e_j^T (output)
i_b = ceil(CasoBs / 3);
j_c = mod(CasoBs - 1, 3) + 1;
BC  = zeros(n);
BC(i_b, j_c) = 1;                                   % B*C

% Mismatch of every node: node 1 nominal, single mismatch on node 2 or 3
param_names = {'da2', 'db2', 'dc2', 'da3', 'db3', 'dc3'};
par_letter  = {'a', 'b', 'c'};
q_idx   = mod(mismatch_idx - 1, 3) + 1;             % a, b or c
node_mm = 2 + (mismatch_idx > 3);                   % node 2 or 3
DP      = zeros(3, 3);                              % row k = mismatch of node k
DP(node_mm, q_idx) = mismatch_ratio * OscPar(q_idx);

% =========================================================================
% STEP 3: SENSITIVITY PROFILE (prediction)
% =========================================================================
if isempty(profile_file)
    cand = dir(fullfile('mat', sprintf('SensitivityNorm_vs_ks_%s_%s_B%d_N*_T*.mat', FMtype, topology, CasoBs)));
    if isempty(cand)
        error('No sensitivity profile for %s, %s, CasoBs = %d in mat/ (run SensitivityNorm_vs_ks.m first).', ...
            FMtype, topology, CasoBs);
    end
    [~, k_new] = max([cand.datenum]);
    profile_file = fullfile(cand(k_new).folder, cand(k_new).name);
end
S = load(profile_file);
fprintf('Sensitivity profile: %s\n', profile_file);
if ~isfield(S, 'ic_mode') || ~strcmpi(S.ic_mode, 'manifold')
    warning('The profile file was not computed with ic_mode = ''manifold''.');
end
ks_prof  = S.ks_sweep(:);
num_pairs = 3;
cols     = (mismatch_idx - 1) * num_pairs + (1:num_pairs);   % e12, e13, e23 of q_j
if ~isfield(S, 'xref')
    error('The profile file has no xref field (run ShiftSensitivity_xref.m first).');
end
xref     = S.xref;                                 % attractor scale x_ref
pred     = xref * 10.^(S.norm_dB(:, cols) / 20) * abs(mismatch_ratio);        % predicted RMS(e_ij)
pred_lo  = xref * 10.^((S.norm_dB(:, cols) - 2*S.norm_unc(:, cols)) / 20) * abs(mismatch_ratio);
pred_hi  = xref * 10.^((S.norm_dB(:, cols) + 2*S.norm_unc(:, cols)) / 20) * abs(mismatch_ratio);
no_infl  = S.no_influence(cols);                    % zero sensitivity vectors

T_end       = S.T_end;
t_transient = S.t_transient;
dt_sample   = S.dt_sample;
t_samp      = t_transient + (1:floor((T_end - t_transient) / dt_sample)) * dt_sample;
X1_pool     = S.X1_pool;
n_runs      = min(n_runs, size(X1_pool, 1));

if isempty(ks_measure)
    ks_measure = ks_prof.';
end
N_meas = numel(ks_measure);

% =========================================================================
% STEP 4: MEASURED RMS OF THE SYNCHRONIZATION ERRORS (full network)
% =========================================================================
P_all = repmat(OscPar, 3, 1) + DP;                  % parameters of nodes 1, 2, 3
rms_runs = NaN(N_meas, n_runs, num_pairs);
ode_opts = odeset('RelTol', ode_RelTol, 'AbsTol', ode_AbsTol, ...
                  'Events', @(t, X) blowupEvent(X, err_bound));
fprintf('Measured RMS: %s, %s, b_%d -> c_%d, %s = %.3g (r = %g), %d ks x %d runs, window [%g, %g]\n', ...
    FMtype, topology, i_b, j_c, param_names{mismatch_idx}, DP(node_mm, q_idx), mismatch_ratio, ...
    N_meas, n_runs, t_transient, T_end);
t_start = tic;
for k = 1:N_meas
    ks_val = ks_measure(k);
    f_net  = @(t, X) networkRHS(X, Fnode, P_all, edges, BC, ks_val, n);
    for r = 1:n_runs
        x1_0 = X1_pool(r, :).';
        X0   = [x1_0; x1_0; x1_0];                  % Start on the synchronization manifold
        [t_out, X_out, te] = ode45(f_net, [0, t_samp], X0, ode_opts);
        if ~isempty(te) || numel(t_out) < numel(t_samp) + 1
            continue;                               % Blow-up: run excluded (NaN)
        end
        Xs = X_out(2:end, :);                       % Samples in the RMS window
        e  = {Xs(:, 1:3) - Xs(:, 4:6), Xs(:, 1:3) - Xs(:, 7:9), Xs(:, 4:6) - Xs(:, 7:9)};
        for p = 1:num_pairs
            rms_runs(k, r, p) = sqrt(sum(mean(e{p}.^2, 1)));   % vector RMS
        end
    end
    fprintf('  ks = %.4f (%d/%d)  median RMS(e12, e13, e23) = %.4g, %.4g, %.4g  (%.1f min)\n', ks_val, k, N_meas, ...
        median(rms_runs(k, :, 1), 'omitnan'), median(rms_runs(k, :, 2), 'omitnan'), ...
        median(rms_runs(k, :, 3), 'omitnan'), toc(t_start)/60);
end
meas = squeeze(median(rms_runs, 2, 'omitnan'));     % N_meas x 3
if N_meas == 1, meas = meas(:).'; end

% Agreement at the ks where both exist (ratio measured / predicted)
[ks_common, i_m, i_p] = intersect(round(ks_measure(:), 10), round(ks_prof, 10));
ratio = meas(i_m, :) ./ pred(i_p, :);
fprintf('\nMeasured / predicted at the common ks values:\n');
fprintf('%10s %12s %12s %12s\n', 'ks', 'e12', 'e13', 'e23');
for k = 1:numel(ks_common)
    fprintf('%10.4f %12.3f %12.3f %12.3f\n', ks_common(k), ratio(k, 1), ratio(k, 2), ratio(k, 3));
end

% =========================================================================
% STEP 5: FIGURE - measured vs predicted, one panel per error vector
% =========================================================================
col_meas = [0 0.447 0.741];
col_pred = [0.850 0.325 0.098];
err_tex  = {'12', '13', '23'};
fig = figure('Name', 'Measured vs predicted RMS error', 'Position', fig_size);
tl  = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for p = 1:num_pairs
    ax = nexttile(tl);
    hold(ax, 'on');
    h = gobjects(1, 2);
    if ~no_infl(p)
        ok = isfinite(pred(:, p));
        if show_band && any(ok)
            fill(ax, [ks_prof(ok); flipud(ks_prof(ok))], [pred_lo(ok, p); flipud(pred_hi(ok, p))], ...
                col_pred, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        h(2) = plot(ax, ks_prof, pred(:, p), '--s', 'Color', col_pred, 'MarkerFaceColor', col_pred, ...
            'MarkerSize', 5, 'LineWidth', 2);
    end
    h(1) = plot(ax, ks_measure, meas(:, p), '-o', 'Color', col_meas, 'MarkerFaceColor', col_meas, ...
        'MarkerSize', 4, 'LineWidth', 2);
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontSize', font_ticks, 'TickLabelInterpreter', 'latex');
    xlim(ax, [0, ceil(max([ks_measure(:); ks_prof]))]);
    if y_log
        set(ax, 'YScale', 'log');                   % Whole decades covering every measured and predicted value
        y_all = [meas(:, p); pred(:, p)];
        y_all = y_all(isfinite(y_all) & y_all > 0);
        if ~isempty(y_all)
            ylim(ax, 10.^[floor(log10(min(y_all))), ceil(log10(max(y_all)))]);
        end
    elseif ~isempty(y_max)
        ylim(ax, [0, y_max]);
    else                                            % No measured curve leaves the panel
        y_top = max(meas(isfinite(meas(:, p)), p));
        if ~isempty(y_top) && y_top > 0
            ylim(ax, [0, ceil(1.1 * y_top / 0.5) * 0.5]);
        end
    end
    ylabel(ax, sprintf('$\\mathrm{RMS}(\\mathbf{e}_{%s})$', err_tex{p}), 'Interpreter', 'latex', 'FontSize', font_label);
    text(ax, 0.97, 0.90, sprintf('(%s)', char('a' + p - 1)), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'Interpreter', 'latex', 'FontSize', font_ticks);
    if no_infl(p)
        legend(ax, h(1), {'Measured (no influence predicted)'}, 'Interpreter', 'latex', ...
            'FontSize', font_legend, 'Location', 'northeast');
    else
        legend(ax, h, {'Measured', 'Predicted'}, 'Interpreter', 'latex', ...
            'FontSize', font_legend, 'Location', 'northeast');
    end
end
xlabel(tl, '$k_s$', 'Interpreter', 'latex', 'FontSize', font_label);
if show_title
    title(tl, sprintf('%s - %s, $\\mathbf{b}_%d \\to \\mathbf{c}_%d$: $\\delta %s_%d = %g\\,%s$', osc_tex{OsID}, ...
        motif_names{topology_idx}, i_b, j_c, par_letter{q_idx}, node_mm, mismatch_ratio, par_letter{q_idx}), ...
        'Interpreter', 'latex', 'FontSize', font_title);
end
all_ax = findall(fig, 'Type', 'axes');
for kk = 1:numel(all_ax)
    all_ax(kk).Toolbar.Visible = 'off';
end

% =========================================================================
% STEP 6: SAVE (eps/ and mat/)
% =========================================================================
if ~exist('eps', 'dir'), mkdir('eps'); end
if ~exist('mat', 'dir'), mkdir('mat'); end
file_tag = sprintf('%s_%s_B%d_%s_r%g', FMtype, topology, CasoBs, param_names{mismatch_idx}, mismatch_ratio);
exportgraphics(fig, fullfile('eps', sprintf('ValidationRMS_vs_ks_%s.eps', file_tag)), 'ContentType', 'vector');
save(fullfile('mat', sprintf('ValidationRMS_vs_ks_%s.mat', file_tag)), ...
    'ks_measure', 'rms_runs', 'meas', 'ks_prof', 'pred', 'pred_lo', 'pred_hi', 'no_infl', ...
    'ks_common', 'ratio', 'mismatch_idx', 'mismatch_ratio', 'DP', 'OscPar', 'OsID', 'topology_idx', ...
    'CasoBs', 'n_runs', 'T_end', 't_transient', 'dt_sample', 'ode_RelTol', 'ode_AbsTol', 'profile_file', 'xref');
fprintf('Saved: eps/ValidationRMS_vs_ks_%s.eps and mat/ValidationRMS_vs_ks_%s.mat\n', file_tag, file_tag);

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function dX = networkRHS(X, Fnode, P_all, edges, BC, ks, n)
    % dx_k/dt = F(x_k; q + dq_k) + ks * sum_{s -> k} B*C*(x_s - x_k)
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

function [value, isterminal, direction] = blowupEvent(X, err_bound)
    if any(~isfinite(X))
        value = -1;
    else
        value = err_bound - max(abs(X));
    end
    isterminal = 1;
    direction  = 0;
end

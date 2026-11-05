% =========================================================================
% SCRIPT: ModalSensitivity_vs_ks.m
% DESCRIPTION: Eigenmode (master stability) prediction of the error
%              sensitivities of the three-node motifs, compared with the
%              direct node-space computation of SensitivityNorm_vs_ks*.m.
%              On the synchronization manifold, with Laplacian
%              L = V diag(mu) V^-1 (Appendix A.1), a mismatch dq_j at node m
%              drives each transverse mode k with the master-slave
%              sensitivity s(t; kappa), kappa = mu_k ks:
%                  ds/dt = (DF(x1) - kappa b c^T) s + J_q(x1),  s(0) = 0
%              and the error sensitivity of e_ij is
%                  s_eij(t) = sum_k (V_ik - V_jk) (V^-1)_km s(t; mu_k ks).
%              The topology enters only through mu_k and the coefficients
%              (V_ik - V_jk)(V^-1)_km; the ks dependence only through the
%              master-slave curve. The rms, the median over the runs and the
%              vector norm follow Eqs. (eq_sensitivity_index_vector) and
%              (eq_Sq_bar), with the same node-1 pool (base_seed), window and
%              solver settings as SensitivityNorm_vs_ks.m, so the prediction
%              and the direct values differ only by the method.
%              Motifs with complex eigenvalues (unidirectional ring chain),
%              a non-diagonalizable Laplacian (linear chain, bidirectional
%              mixed ring) or no stable domain (linear competitive) are
%              skipped.
%              Input : mat/SensitivityNorm_vs_ks_<Osc>_AllTopologies_N<N>_T<T_end>.mat
%              Output: mat/ModalSensitivity_<Osc>_B<C>_N<N>_T<T_end>.mat
%                      eps/ModalSensitivity_<Osc>_B<C>_N<N>_T<T_end>.eps
%                      tex/ModalSensitivity_<Osc>_B<C>_N<N>_T<T_end>.tex
% =========================================================================
clear; clc; close all;

% =========================================================================
% STEP 1: CONFIGURATION
% =========================================================================
OsID         = 1;               % 1=Rossler, 2=Lorenz
CasoBs       = 1;               % 1..9 (b_i -> c_j, SISObc)
N_sweep      = 25;              % Same N_sweep as the computation
time_profile = 'medium';        % Same time_profile as the computation
topology_run = [];              % [] = every diagonalizable motif with a result
n_runs       = 20;              % Realizations (node-1 pool points 1..n_runs)

dt_sample        = 0.05;        % Same settings as SensitivityNorm_vs_ks.m
ode_RelTol       = 1e-6;
ode_AbsTol       = 1e-8;
base_seed        = 42;
pool_T_transient = 500;
pool_spacing     = 50;
cond_max         = 1e8;         % cond(V) above this: Laplacian treated as not diagonalizable
zero_coef_tol    = 1e-9;        % modal coefficients below this are zero

mat_folder = 'mat';
eps_folder = 'eps';
tex_folder = 'tex';
show_title = true;
font_title = 14; font_label = 14; font_ticks = 12; font_legend = 11;

OscType    = {'Rossler', 'Lorenz', 'Chen', 'HR', 'Duffing', 'vanderPol'};
topo_short = {'Lin. chain', 'Lin. aux.', 'Lin. relay', 'Lin. mix. mast.', 'Lin. comp.', ...
    'Lin. mix. slav.', 'Ring bidir. chain', 'Ring bidir. dbl. mix.', 'Ring bidir. mix.', ...
    'Ring unidir. chain', 'Ring bidir. aux.', 'Ring unidir. aux.', 'Ring bidir. comp.'};
% Edge lists [i j] = node i drives node j (same as MotifLLE_vs_ks_AllTopologies.m)
edges_all = {[1 2; 2 3], [1 2; 1 3], [1 2; 2 1; 2 3; 3 2], [1 2; 2 3; 3 2], [1 2; 3 2], ...
    [2 1; 2 3; 3 2], [1 2; 2 1; 2 3; 3 2; 3 1; 1 3], [1 2; 2 3; 3 2; 3 1; 1 3], ...
    [1 2; 2 3; 3 2; 3 1], [1 2; 2 3; 3 1], [1 2; 1 3; 2 3; 3 2], [1 2; 1 3; 2 3], ...
    [1 2; 2 1; 2 3; 1 3]};
switch lower(time_profile)
    case 'short',  t_transient = 200;  T_end = 1000;
    case 'medium', t_transient = 500;  T_end = 2500;
    case 'long',   t_transient = 1000; T_end = 5000;
    otherwise, error('time_profile must be ''short'', ''medium'' or ''long''.');
end
pairs = [1 2; 1 3; 2 3];        % e12, e13, e23 (v = (j-1)*3 + p, as in SensitivityNorm_vs_ks.m)

% =========================================================================
% STEP 2: OSCILLATOR, MASTER-SLAVE SENSITIVITY FUNCTIONS
% =========================================================================
syms a b c real
[F, OscPar, FMtype] = OscillatorType(OsID);
syms x(t) y(t) z(t)
xs = sym('xs', [3 1], 'real');
Fs = subs(F, [x(t); y(t); z(t)], xs);
n  = 3;
DFs = jacobian(Fs, xs);                                % n x n
Jqs = jacobian(Fs, [a b c]);                           % n x 3
nom = OscPar(1:3);
f_x  = matlabFunction(subs(Fs,  [a b c], nom), 'Vars', {xs});
f_DF = matlabFunction(subs(DFs, [a b c], nom), 'Vars', {xs});
f_Jq = matlabFunction(subs(Jqs, [a b c], nom), 'Vars', {xs});
[bs, cs] = SISObc(CasoBs, n);
BC = bs * cs;                                          % b_p c_q^T
xref = AttractorScale(OsID);
p_scale = nom;                                         % sbar = s q_j / x_ref

% Node-1 pool: identical to SensitivityNorm_vs_ks.m
rng(base_seed, 'twister');
[~, Xp] = ode45(@(t, x) f_x(x), [0, pool_T_transient + (1:n_runs) * pool_spacing], ...
                rand(n, 1), odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
X1_pool = Xp(2:end, :);
t_samp  = t_transient + (1:floor((T_end - t_transient) / dt_sample)) * dt_sample;

% =========================================================================
% STEP 3: DIRECT RESULTS
% =========================================================================
filename_mat = fullfile(mat_folder, sprintf('SensitivityNorm_vs_ks_%s_AllTopologies_N%d_T%d.mat', FMtype, N_sweep, T_end));
if ~isfile(filename_mat)
    error('%s not found: run SensitivityNorm_Consolidate.m first.', filename_mat);
end
U = load(filename_mat, 'SN', 'TS_tab');
has_result = any(isfinite(U.TS_tab(:, CasoBs, :)), 3);
if isempty(topology_run), topology_run = find(has_result).'; end

% =========================================================================
% STEP 4: MODAL PREDICTION PER MOTIF
% =========================================================================
res = struct([]);
for ti = topology_run
    if ~has_result(ti), continue; end
    E = edges_all{ti};
    A = zeros(3);
    for e = 1:size(E, 1), A(E(e, 2), E(e, 1)) = 1; end  % row = driven node
    L = diag(sum(A, 2)) - A;
    [V, D] = eig(L);
    mu = diag(D);
    if any(abs(imag(mu)) > 1e-9)
        fprintf('%s: complex eigenvalues, skipped\n', topo_short{ti}); continue;
    end
    mu = real(mu); V = real(V);
    if cond(V) > cond_max
        fprintf('%s: Laplacian not diagonalizable, skipped\n', topo_short{ti}); continue;
    end
    if nnz(abs(mu) < 1e-9) > 1
        fprintf('%s: repeated zero eigenvalue, skipped\n', topo_short{ti}); continue;
    end
    Vi = inv(V);
    kmod = find(abs(mu) > 1e-9).';                     % the two transverse modes
    % Modal coefficients c(p, m, k): error pair p, mismatch node m (2,3), mode k
    coef = zeros(3, 3, 2);
    for p = 1:3
        for m = 2:3
            for kk = 1:2
                k = kmod(kk);
                coef(p, m, kk) = (V(pairs(p,1), k) - V(pairs(p,2), k)) * Vi(k, m);
            end
        end
    end
    coef(abs(coef) < zero_coef_tol) = 0;

    R  = U.SN(ti, CasoBs);
    ks = R.ks_sweep(:);
    Nk = numel(ks);
    kap = [mu(kmod(1)) * ks; mu(kmod(2)) * ks];       % 2 Nk couplings, mode-major
    Nkap = numel(kap);

    % Per run: rms of each component of sbar_eij, for every ks, vector v, component
    RMS = NaN(n_runs, Nk, 18, n);
    fprintf('%s: mu = [%.4g %.4g], %d couplings, %d runs ', topo_short{ti}, mu(kmod), Nkap, n_runs);
    for r = 1:n_runs
        Z0 = [X1_pool(r, :).'; zeros(n * 3 * Nkap, 1)];
        [~, Z] = ode45(@(t, Z) modalRHS(Z, f_x, f_DF, f_Jq, BC, kap, n), [0, t_samp], Z0, ...
                       odeset('RelTol', ode_RelTol, 'AbsTol', ode_AbsTol));
        Z = Z(2:end, :);                                % samples in the window
        Sk = reshape(Z(:, n+1:end), [], n, 3, Nkap);    % time x component x parameter x kappa
        for jj = 1:3                                    % parameter a, b, c
            for m = 2:3                                 % node of the mismatch
                j = (m - 2) * 3 + jj;                   % q_j index 1..6
                for p = 1:3
                    cf = squeeze(coef(p, m, :));
                    if all(cf == 0), continue; end
                    v = (j - 1) * 3 + p;
                    for kk = 1:Nk
                        s_e = cf(1) * Sk(:, :, jj, kk) + cf(2) * Sk(:, :, jj, Nk + kk);
                        RMS(r, kk, v, :) = sqrt(mean(s_e.^2, 1)) * p_scale(jj) / xref;
                    end
                end
            end
        end
        fprintf('.');
    end
    fprintf('\n');
    med = squeeze(median(RMS, 1));                      % Nk x 18 x n (linear medians)
    norm_pred = 20 * log10(sqrt(sum(med.^2, 3)));       % Nk x 18, NaN where structurally zero
    % Sbar_q(e) curve, Eq. (eq_Sq_bar): zero vectors count as 0
    S_lin = 10.^(norm_pred / 20); S_lin(~isfinite(S_lin)) = 0;
    TS_pred = NaN(Nk, 6);
    for j = 1:6
        TS_pred(:, j) = 20 * log10(sqrt(sum(S_lin(:, (j-1)*3 + (1:3)).^2, 2) / 3));
    end
    norm_dir = R.norm_dB;
    TS_dir   = R.TS_curve_dB;
    d_vec = norm_pred - norm_dir;   d_vec(~isfinite(d_vec)) = NaN;
    d_TS  = TS_pred - TS_dir;       d_TS(~isfinite(d_TS)) = NaN;

    k = numel(res) + 1;
    res(k).ti = ti;              res(k).name = topo_short{ti};
    res(k).mu = mu(kmod).';      res(k).coef = coef;
    res(k).ks = ks;              res(k).norm_pred = norm_pred;  res(k).norm_dir = norm_dir;
    res(k).TS_pred = TS_pred;    res(k).TS_dir = TS_dir;
    res(k).avg_pred = mean(TS_pred, 1, 'omitnan');
    res(k).avg_dir  = mean(TS_dir, 1, 'omitnan');
    res(k).rms_vec  = sqrt(mean(d_vec(:).^2, 'omitnan'));
    res(k).max_vec  = max(abs(d_vec(:)), [], 'omitnan');
    res(k).rms_TS   = sqrt(mean(d_TS(:).^2, 'omitnan'));
    res(k).max_TS   = max(abs(d_TS(:)), [], 'omitnan');
    res(k).max_avg  = max(abs(res(k).avg_pred - res(k).avg_dir));
    fprintf('   vectors: rms %.2f, max %.2f dB | Sbar_q: rms %.2f, max %.2f dB | <Sbar_q>: max %.2f dB\n', ...
        res(k).rms_vec, res(k).max_vec, res(k).rms_TS, res(k).max_TS, res(k).max_avg);
end

filename_out = fullfile(mat_folder, sprintf('ModalSensitivity_%s_B%d_N%d_T%d.mat', FMtype, CasoBs, N_sweep, T_end));
save(filename_out, 'res', 'CasoBs', 'n_runs', 'X1_pool', 'xref', 'filename_mat', 'FMtype');
fprintf('Saved: %s\n', filename_out);

% =========================================================================
% STEP 5: FIGURE: Sbar_q(e) vs ks, direct (lines) and modal (markers)
% =========================================================================
if ~exist(eps_folder, 'dir'), mkdir(eps_folder); end
nm = numel(res);
fig = figure('Units', 'pixels', 'Position', [100, 100, 900, 260 * ceil(nm / 3)]);
tl = tiledlayout(fig, ceil(nm / 3), 3, 'TileSpacing', 'compact', 'Padding', 'compact');
cols = lines(6);
par_tex = {'\delta a_2', '\delta b_2', '\delta c_2', '\delta a_3', '\delta b_3', '\delta c_3'};
for k = 1:nm
    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'on');
    hl = gobjects(1, 6);
    for j = 1:6
        hl(j) = plot(ax, res(k).ks, res(k).TS_dir(:, j), '-', 'Color', cols(j, :), 'LineWidth', 1.2);
        plot(ax, res(k).ks, res(k).TS_pred(:, j), 'o', 'Color', cols(j, :), 'MarkerSize', 4);
    end
    title(ax, sprintf('%s ($\\mu = %.3g, %.3g$)', res(k).name, res(k).mu), 'Interpreter', 'latex', 'FontSize', font_title);
    xlabel(ax, '$k_s$', 'Interpreter', 'latex', 'FontSize', font_label);
    ylabel(ax, '$\bar{S}_q(\mathbf{e})$ (dB)', 'Interpreter', 'latex', 'FontSize', font_label);
    ax.TickLabelInterpreter = 'latex'; ax.FontSize = font_ticks;
end
lg = legend(ax, hl, strcat('$', par_tex, '$'), 'Interpreter', 'latex', 'FontSize', font_legend, 'NumColumns', 6);
lg.Layout.Tile = 'south';
if show_title
    [p_, q_] = ind2sub([3 3], CasoBs);
    title(tl, sprintf('%s, $\\mathbf{b}_%d \\to \\mathbf{c}_%d$: direct (lines) and eigenmode prediction (markers)', ...
        FMtype, q_, p_), 'Interpreter', 'latex', 'FontSize', font_title);
end
filename_eps = fullfile(eps_folder, sprintf('ModalSensitivity_%s_B%d_N%d_T%d.eps', FMtype, CasoBs, N_sweep, T_end));
exportgraphics(fig, filename_eps, 'ContentType', 'vector');
close(fig);
fprintf('Saved: %s\n', filename_eps);

% =========================================================================
% STEP 6: LATEX TABLE OF THE DIFFERENCES -> tex/
% =========================================================================
if ~exist(tex_folder, 'dir'), mkdir(tex_folder); end
filename_tex = fullfile(tex_folder, sprintf('ModalSensitivity_%s_B%d_N%d_T%d.tex', FMtype, CasoBs, N_sweep, T_end));
fid = fopen(filename_tex, 'w');
fprintf(fid, '%% Generated by ModalSensitivity_vs_ks.m from %s\n', filename_mat);
fprintf(fid, '\\begin{tabular}{lcccc}\n\\toprule\n');
fprintf(fid, 'Motif & $\\mu_2, \\mu_3$ & $S_q(\\mathbf{e}_{ij})$ rms / max (dB) & $\\bar{S}_q(\\mathbf{e})$ rms / max (dB) & $\\langle \\bar{S}_q(\\mathbf{e})\\rangle$ max (dB) \\\\\n\\midrule\n');
for k = 1:nm
    fprintf(fid, '%s & $%.3g,\\ %.3g$ & $%.2f$ / $%.2f$ & $%.2f$ / $%.2f$ & $%.2f$ \\\\\n', res(k).name, res(k).mu, ...
        res(k).rms_vec, res(k).max_vec, res(k).rms_TS, res(k).max_TS, res(k).max_avg);
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);
fprintf('Saved: %s\n', filename_tex);

% =========================================================================
% LOCAL FUNCTION: master state + master-slave sensitivities at every kappa
%   Z = [x1; vec(S_1); ...; vec(S_K)], S_k = n x 3 (columns a, b, c)
% =========================================================================
function dZ = modalRHS(Z, f_x, f_DF, f_Jq, BC, kap, n)
    x1 = Z(1:n);
    S  = reshape(Z(n+1:end), n, 3, []);                 % n x 3 x K
    DF = f_DF(x1);
    Jq = f_Jq(x1);
    dS = pagemtimes(DF, S) - pagemtimes(BC, S) .* reshape(kap, 1, 1, []) + Jq;
    dZ = [f_x(x1); dS(:)];
end

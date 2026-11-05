% =========================================================================
% SCRIPT: SensitivityNorm_vs_ks.m
% DESCRIPTION: Self-contained. RMS sensitivity vector index ||s_qj(e_ij)||
%              (dB) of the motif error system against the coupling
%              strength ks, for ONE topology and ONE coupling scheme, with
%              ks swept inside a synchronization interval (lambda_max < 0)
%              read from mat/LLE_vs_ks_<Osc>_AllTopologies.mat. The
%              interval is trimmed to the ks where the transverse
%              transient has died out before the RMS window starts,
%              |lambda_max| * t_transient >= n_decay, and (use_rms_file)
%              to rms(e) < rms_file_tol from mat/RMS_vs_ks_<Osc>_<topology>_T<T_end>.mat,
%              before any run. Optionally
%              (use_internal_rms_check) every ks whose synchronization
%              error in the same runs, median over runs of rms(e) in the
%              RMS window, exceeds rms_e_tol (e.g. bubbling near a
%              threshold) is flagged and left out of the figures.
%              Part 1 builds the sensitivity matrix equation of the
%              mismatched error system (only the code needed for it) and
%              generates odeMotifSensitivityMismatch.m, with ks as a
%              runtime parameter (one file for the whole sweep).
%              Part 2, for each of the N_sweep values of ks: runs in batches
%              until the same convergence test as SensitivityIndexMonteCarlo.m is met
%              (start on the synchronization manifold: node 1 from an
%              attractor pool, e(0) = 0, S(0) = 0, sensitivity of node 1
%              held at 0; ALL mismatches = 0), ode45, RMS of
%              sbar_ij = s_ij * p_j / x_ref for t >= t_transient; per component
%              the median over runs (dB), per error vector the norm
%              (Eq. 23, on linear values) in dB and its uncertainty
%              1.2533 * std(per-run norm in dB) / sqrt(n). The same pool
%              and seeds are used at every ks (common random numbers), so
%              the curves are smooth in ks.
%              Figures (eps/):
%                1) ||s_qj(e)|| (dB) vs ks: 2x3 tiles, one per mismatch,
%                   three curves e12, e13, e23 (+/- 2 x uncertainty);
%                2) Robustness summary: (a) lambda_max(ks) with the
%                   swept interval and the most robust ks; (b) worst-case
%                   and mean sensitivity vs ks; (c) sensitivity against
%                   |lambda_max| (log-log) with the linear-response slope
%                   ||s|| ~ 1/|lambda_max|.
%              Results saved to mat/.
% =========================================================================
clear; clc; close all;

% =========================================================================
% CONFIGURATION
% =========================================================================
OsID         = 1;               % 1=Rossler, 2=Lorenz
topology_idx = 1;               % 1..13 (list below)
CasoBs       = 5;               % 1..9 (b_i -> c_j, SISObc)

% --- ks sweep inside a synchronization interval of the LLE file ---
N_sweep       = 25;             % Number of ks values
interval_pick = 0;              % 0 = longest interval; k = k-th interval (in ks order)
% Edge trimming (lambda criterion): keep only ks whose transverse
% transient has died out before the RMS window starts,
%   |lambda_max| * t_transient >= n_decay   (error decayed by e^n_decay)
n_decay       = 5;              % 5 -> |lambda_max| >= 5 / t_transient
% Double check with the synchronization error of the SAME runs (same
% window and initial conditions as the index): a ks is flagged, and left
% out of the figures, when the median over runs of rms(e) for
% t >= t_transient exceeds rms_e_tol, e = [e12; e13; e23]
use_internal_rms_check = true;
rms_e_tol              = 1e-3;
% Pre-trim with the RMS synchronization error already computed in
% mat/RMS_vs_ks_<Osc>_<topology>_T<T_end>.mat (Motif_RMS_vs_ks_AllTopologies.m run with the
% SAME time_profile, so both use the same window):
% keep only ks with rms(e) < rms_file_tol before any run is spent. That file
% may use a shorter horizon than this script: slowly converging ks near an
% edge can fail it although the lambda criterion accepts them.
use_rms_file = true;
rms_file_tol = [];              % [] = sync_tol stored in the RMS file

% --- Runs per ks ---
% Runs per ks: same convergence test as SensitivityIndexMonteCarlo.m,
% applied at every ks (median index of the 54 coefficients)
% Safety cap chosen from the oscillator (only reached if a ks never converges):
%              Rossler  Lorenz  Chen   HR    Duffing  vanderPol
max_runs_table = [ 100,  200,    200,   200,  200,     200 ];
max_runs     = max_runs_table(OsID);
batch_size   = 10;              % Runs between convergence checks
min_runs     = 10;              % Valid runs before the test is applied
conv_window  = 10;              % Compare with the index conv_window runs earlier
conv_tol_dB  = 1;               % max |change of the index| over the window (dB)
conv_sem_dB  = 1;               % max uncertainty of the index (dB)
% true: no runs; the profile mat/SensitivityNorm_vs_ks_<tag>.mat (with xref)
% is loaded and only the figures are redrawn (the .mat file is not rewritten)
replot_only  = false;
% Time horizon profile (sets t_transient and T_end; the lambda threshold
% lambda_min = n_decay / t_transient follows automatically):
%   'short'  : t_transient =  200, T_end = 1000  -> lambda_min = 0.025  (1x cost)
%   'medium' : t_transient =  500, T_end = 2500  -> lambda_min = 0.010  (2.5x)
%   'long'   : t_transient = 1000, T_end = 5000  -> lambda_min = 0.005  (5x)
time_profile = 'medium';
switch lower(time_profile)
    case 'short',  t_transient = 200;  T_end = 1000;
    case 'medium', t_transient = 500;  T_end = 2500;
    case 'long',   t_transient = 1000; T_end = 5000;
    otherwise, error('time_profile must be ''short'', ''medium'' or ''long''.');
end
dt_sample   = 0.05;             % Sampling step of the RMS
ode_RelTol  = 1e-6;
ode_AbsTol  = 1e-8;
base_seed   = 42;               % Run r uses rng(base_seed + r) (same at every ks)
pool_T_transient = 500;         % Node-1 pool: discarded transient
pool_spacing     = 50;          % Node-1 pool: time between stored points
err_bound    = 1e3;             % Blow-up bound on |X1| and |e_ij| (run excluded)
zero_rel_tol = 1e-6;            % Zero coefficient: mean RMS <= zero_rel_tol * max at that ks

show_errorbars = true;          % +/- 2 x uncertainty on the norm curves (figure 1)
marker_step    = 3;             % Figure 1: a marker every marker_step points, staggered between curves
font_title    = 20;
font_subtitle = 20;
font_xlabel   = 20;
font_ylabel   = 20;
font_legend   = 12;
font_ticks    = 15;
show_title    = true;           % Figure 1: false -> no title (the paper caption describes the case)
show_panel_letters = false;     % Figure 1: (a)-(f) inside the panels; false when the figure is a subfigure of the paper
fig1_size     = [100 100 1000 650];   % Figure 1 size (pixels): room for labels and legends

% #########################################################################
% PART 1: SENSITIVITY MATRIX EQUATION -> odeMotifSensitivityMismatch.m
% #########################################################################
% =========================================================================
% STEP 1: INITIALIZATION & OSCILLATOR SETUP
% =========================================================================
num_oscillators = 3;
syms a b c ks t real
[F, OscPar, FMtype] = OscillatorType(OsID);
syms x(t) y(t) z(t)
X = [x(t); y(t); z(t)];
n = length(X);          % State dimension per node (3)

% --- Parameter mismatches of nodes 2 and 3 (node 1 nominal) ---
syms da2 db2 dc2
syms da3 db3 dc3
P   = [a b c];
DP2 = [da2 db2 dc2];
DP3 = [da3 db3 dc3];
% DP_all{k} is the mismatch applied to oscillator k. Node 1 = zero
% mismatch (nominal), so the same subs() call can run for every k
% uniformly without an if/else in the loop.
DP_all = {sym(zeros(1, numel(P))), DP2, DP3};
Par = [P DP2 DP3 ks];
total_states = n * num_oscillators;
num_pairs  = (num_oscillators * (num_oscillators - 1)) / 2; % Unique error pairs
dim_Y      = n * num_pairs;
% =========================================================================
% STEP 2: TOPOLOGY SELECTION & COUPLING SETUP
% =========================================================================
edges_chain                   = [1 2; 2 3];
edges_auxiliar                = [1 2; 1 3];
edges_relay                   = [1 2; 2 1; 2 3; 3 2];
edges_mixed_master            = [1 2; 2 3; 3 2];
edges_competitive             = [1 2; 3 2];
edges_mixed_slave             = [2 1; 2 3; 3 2];
edges_ring_bidir_chain        = [1 2; 2 1; 2 3; 3 2; 3 1; 1 3];
edges_ring_bidir_double_mixed = [1 2; 2 3; 3 2; 3 1; 1 3];
edges_ring_bidir_mixed        = [1 2; 2 3; 3 2; 3 1];
edges_ring_unidir_chain       = [1 2; 2 3; 3 1];
edges_ring_bidir_auxiliar     = [1 2; 1 3; 2 3; 3 2];
edges_ring_unidir_auxiliar    = [1 2; 1 3; 2 3];
edges_ring_bidir_competitive  = [1 2; 2 1; 2 3; 1 3];
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
edges_input = topology_list{topology_idx};
switch char(edges_input)
case 'chain',                       edges = edges_chain;
case 'auxiliar',                    edges = edges_auxiliar;
case 'relay',                       edges = edges_relay;
case 'mixed_master',                edges = edges_mixed_master;
case 'competitive',                 edges = edges_competitive;
case 'mixed_slave',                 edges = edges_mixed_slave;
case 'ring_bidir_chain',            edges = edges_ring_bidir_chain;
case 'ring_bidir_double_mixed',     edges = edges_ring_bidir_double_mixed;
case 'ring_bidir_mixed',            edges = edges_ring_bidir_mixed;
case 'ring_unidir_chain',           edges = edges_ring_unidir_chain;
case 'ring_bidir_auxiliar',         edges = edges_ring_bidir_auxiliar;
case 'ring_unidir_auxiliar',        edges = edges_ring_unidir_auxiliar;
case 'ring_bidir_competitive',      edges = edges_ring_bidir_competitive;
end
num_edges = size(edges, 1);
% ks is a runtime parameter of the generated ODE file (no ks_val here)
[bs, cs] = SISObc(CasoBs, n);         % Inner coupling vectors
% =========================================================================
% STEP 3: AUTOMATED STATE & VECTOR FIELD GENERATION
% =========================================================================
XAll = [];
FX_All = [];
for k = 1:num_oscillators
    eval(sprintf('syms x%d(t) y%d(t) z%d(t)', k, k, k));
    eval(sprintf('X%d = [x%d(t); y%d(t); z%d(t)];', k, k, k, k));
    eval(sprintf('FX%d = subs(F, X, X%d);', k, k));
    eval(sprintf('FX%d = subs(FX%d, P, P + DP_all{%d});', k, k, k));
    XAll  = [XAll;  eval(sprintf('X%d', k))];
    FX_All = [FX_All; eval(sprintf('FX%d', k))];
end
Xj  = XAll;             % Names used by the rest of the script
FXj = FX_All;
% =========================================================================
% STEP 4: GLOBAL BLOCK COUPLING MATRIX CONSTRUCTION
% =========================================================================
Mc = sym(zeros(total_states, total_states));
idx = cell(num_oscillators, 1);
for k = 1:num_oscillators
    idx{k} = (k-1)*n + (1:n);
end
for e = 1:num_edges
    src = edges(e,1);
    tgt = edges(e,2);
    Mc(idx{tgt}, idx{src}) = Mc(idx{tgt}, idx{src}) + ks * (bs * cs);
    Mc(idx{tgt}, idx{tgt}) = Mc(idx{tgt}, idx{tgt}) - ks * (bs * cs);
end
F_coupling = FXj + Mc * Xj;
% =========================================================================
% DYNAMIC ALL-PAIRS ERROR SYMBOLIC GENERATION (ARRAYS ONLY)
% =========================================================================
E_array    = cell(num_pairs, 1);
pair_names = cell(num_pairs, 1);
pair_ij    = zeros(num_pairs, 2);      % (i, j) of each pair
vars_E     = [];
pair_idx   = 1;
for i = 1:num_oscillators
    for j = (i+1):num_oscillators
        pair_str = sprintf('e%d%d', i, j);
        pair_names{pair_idx} = pair_str;
        pair_ij(pair_idx, :) = [i, j];
        eval(sprintf('syms %sx(t) %sy(t) %sz(t)', pair_str, pair_str, pair_str));
        eval(sprintf('%s = [%sx(t); %sy(t); %sz(t)];', pair_str, pair_str, pair_str, pair_str));
        eval(sprintf('E_array{%d} = %s;', pair_idx, pair_str));
        eval(sprintf('vars_E = [vars_E; %s];', pair_str));
        pair_idx = pair_idx + 1;
    end
end
% =========================================================================
% ERROR SYSTEM:  de_ij/dt = [F_i(X_i) - F_j(X_i - e_ij)] + u_i - u_j
%   F_k = F with parameters P + DP_k (node k's own parameters).
%   u_k = ks * sum_{s -> k} B*C*(X_s - X_k)  is the coupling input of node k,
%   and every difference X_s - X_k is written with its own error symbol
%   (+e_sk if s < k, -e_ks if s > k). Each equation keeps e12, e13, e23;
%   the nonlinear part depends on e_ij, the reference node X_i and the
%   mismatches. With mismatch, e = 0 is no longer a solution: the errors
%   stay bounded near zero at best (practical synchronization).
% =========================================================================

% Pair index lookup: Pidx(i,j) = position of e_ij in E_array
% (named Pidx because P holds the nominal parameters [a b c])
Pidx = zeros(num_oscillators);
for p = 1:num_pairs
    Pidx(pair_ij(p,1), pair_ij(p,2)) = p;
end
errDiff = @(s, k) (s < k) * E_array{Pidx(min(s,k), max(s,k))} ...
                - (s > k) * E_array{Pidx(min(s,k), max(s,k))};   % X_s - X_k

% Coupling input of each node in error variables
U = cell(num_oscillators, 1);
for k = 1:num_oscillators
    U{k} = sym(zeros(n, 1));
end
for e = 1:num_edges
    src = edges(e,1);
    tgt = edges(e,2);
    U{tgt} = U{tgt} + ks * (bs * cs) * errDiff(src, tgt);
end

% Error equations
F_matrix_E = [];
for p = 1:num_pairs
    i = pair_ij(p, 1);
    j = pair_ij(p, 2);
    Xi = Xj(idx{i});
    Fj = subs(F, P, P + DP_all{j});                        % node j's vector field
    FXi_minus = subs(Fj, X, Xi - E_array{p});              % F_j(X_i - e_ij)
    dEij = (FXj(idx{i}) - FXi_minus) + U{i} - U{j};        % F_i(X_i) - F_j(X_i - e_ij)
    F_matrix_E = [F_matrix_E; dEij];
end
F_matrix_E = expand(F_matrix_E);
for r = 1:length(F_matrix_E)
    F_matrix_E(r) = collect(F_matrix_E(r), Par);
end

% =========================================================================
% SUBSTITUTION OF e23 IN de13/dt:  e23 = e13 - e12
% =========================================================================
F_matrix_E(4:6) = expand(subs(F_matrix_E(4:6), E_array{3}, E_array{2} - E_array{1}));
for r = 4:6
    F_matrix_E(r) = collect(F_matrix_E(r), Par);
end

% =========================================================================
% AUGMENTED CLOSED SYSTEM [X1; e12; e13; e23] AND ODE FUNCTION
%   Node 1 supplies the reference trajectory; X2 = X1 - e12 in de23/dt
% =========================================================================
X1v      = Xj(idx{1});
F_X1     = FXj(idx{1}) + U{1};                              % dX1/dt, coupling in errors
F_aug    = [F_X1; F_matrix_E];
F_aug    = subs(F_aug, Xj(idx{2}), X1v - E_array{1});       % X2 = X1 - e12
vars_aug = [X1v; vars_E];

% =========================================================================
% INDEPENDENT SYSTEM x = [X1; e12; e13] AND NOMINAL PARAMETERS
%   de23/dt is redundant (e23 = e13 - e12). Keeping it would add a
%   fictitious-node mode to F_x that round-off can excite. Instead:
%   - substitute e23 = e13 - e12 in dX1/dt, de12/dt and de13/dt,
%   - drop de23/dt (rows 10:12),
%   - substitute the nominal a, b, c numerically (OscPar).
% =========================================================================
F_red    = subs(F_aug(1:3*n), E_array{3}, E_array{2} - E_array{1});
F_red    = subs(F_red, P, OscPar(1:3));             % Nominal a, b, c
vars_red = [X1v; E_array{1}; E_array{2}];           % [X1; e12; e13]

% Safety check: no e23 or X2 may remain in the reduced system
if any(has(F_red, [E_array{3}; Xj(idx{2}); Xj(idx{3})]))
    error('Reduced system still depends on e23, X2 or X3.');
end

% Runtime parameters of the generated file: mismatches and ks
Par_run = [DP2 DP3 ks];

% =========================================================================
% SENSITIVITY MATRIX EQUATION:  dS/dt = F_x(t) S(t) + F_q(t)
%   Independent part: x_I = [X1; e12; e13] (9 states), q = [DP2 DP3] (6)
%     dS_I/dt = F_x S_I + F_q,   F_x = D_x f (9 x 9), F_q = D_q f (9 x 6)
%   e23 part (no feedback, so no fictitious mode):
%     de23/dt  = de13/dt  - de12/dt
%     dS_23/dt = dS_13/dt - dS_12/dt
%   Full state x = [X1; e12; e13; e23] (12), S = dx/dq (12 x 6),
%   s_{i,j}(t) = d x_i / d q_j.
% =========================================================================
q   = [DP2 DP3];                        % Mismatch parameters (sensitivity w.r.t. these)
NxI = numel(vars_red);                  % 9 independent states
Nx  = NxI + n;                          % 12 states including e23
Nq  = numel(q);                         % 6 parameters

% Sensitivity matrix S(t): s_i_j(t) = d x_i / d q_j  (LaTeX: s_{i,j})
syms s(t) [Nx Nq]
S   = s(t);                             % 12 x 6: s1_1(t) ... s12_6(t)
S_I = S(1:NxI, :);                      % Rows of X1, e12, e13

F_x   = jacobian(F_red, vars_red);      % State Jacobian      (9 x 9)
F_q   = jacobian(F_red, q);             % Parametric Jacobian (9 x 6)
F_S_I = F_x * S_I + F_q;                % dS_I/dt             (9 x 6)

% e23 and its sensitivity, driven by the independent rows only
i12 = n + (1:n);                        % rows of e12 in x
i13 = 2*n + (1:n);                      % rows of e13 in x
F_e23  = F_red(i13) - F_red(i12);       % de23/dt  = de13/dt  - de12/dt
F_S_23 = F_S_I(i13, :) - F_S_I(i12, :); % dS_23/dt = dS_13/dt - dS_12/dt

F_x_all = [F_red; F_e23];               % 12 x 1
F_S     = [F_S_I; F_S_23];              % 12 x 6
vars_x  = [vars_red; E_array{3}];       % [X1; e12; e13; e23]

% =========================================================================
% ODE FUNCTION: odeMotifSensitivityMismatch.m
%   dZ = odeMotifSensitivityMismatch(t, Z, Parval)
%   Z      = [x (12); vec(S) (72)], x = [X1; e12; e13; e23],
%            S = reshape(Z(13:end), 12, 6)
%            (column-major: s_1 = dx/dq_1 first, then s_2, ...)
%   Parval = [da2 db2 dc2 da3 db3 dc3 ks]   (same order as Par_run)
%   a, b, c are fixed at OscPar; ks and the mismatches are runtime
%   parameters, so they can be swept without regenerating the file.
% =========================================================================
vars_sens = [vars_x; S(:)];
F_sens    = [F_x_all; F_S(:)];
Nz        = numel(vars_sens);

filename_ode = 'odeMotifSensitivityMismatch.m';
fprintf('Generating MATLAB function file: %s (%d states)...\n', filename_ode, Nz);
odeFunction(F_sens, vars_sens, Par_run, 'File', filename_ode, 'Optimize', true);
clear odeMotifSensitivityMismatch;
rehash;

% #########################################################################
% PART 2: ks SWEEP INSIDE A SYNCHRONIZATION INTERVAL
% #########################################################################

% =========================================================================
% STEP A: SYNCHRONIZATION INTERVAL FROM THE LLE FILE
% =========================================================================
lle_file = fullfile('mat', sprintf('LLE_vs_ks_%s_AllTopologies.mat', FMtype));
if ~isfile(lle_file)
    error('%s not found (run MotifLLE_vs_ks_AllTopologies.m / LLE_Consolidate.m first).', lle_file);
end
L = load(lle_file, 'LLE_array', 'ks_vec', 'kf');
ks_lle  = L.ks_vec(:).';
lle_row = squeeze(L.LLE_array(topology_idx, CasoBs, :)).';
if all(isnan(lle_row))
    error('No LLE data for topology %d, CasoBs %d in %s.', topology_idx, CasoBs, lle_file);
end
mask  = lle_row < 0;                                % NaN (blow-up) counts as not synchronized
d     = diff([0, mask, 0]);
i_on  = find(d == 1);
i_off = find(d == -1) - 1;
keep  = i_off > i_on;                               % Drop single-point intervals
i_on  = i_on(keep);
i_off = i_off(keep);
if isempty(i_on)
    error('No synchronization interval (lambda_max < 0) for %s, b_%d -> c_%d.', ...
        motif_names{topology_idx}, find(bs), find(cs));
end
I_all = [ks_lle(i_on).', ks_lle(i_off).'];
fprintf('\nSynchronization intervals (%s, %s, b_%d -> c_%d):\n', FMtype, motif_names{topology_idx}, find(bs), find(cs));
for r = 1:size(I_all, 1)
    fprintf('  %d: [%.3f, %.3f]%s\n', r, I_all(r,1), I_all(r,2), ternaryText(i_off(r) == numel(ks_lle), '  (reaches kf)', ''));
end
if interval_pick == 0
    [~, r_pick] = max(I_all(:,2) - I_all(:,1));
else
    r_pick = interval_pick;
end
fprintf('Chosen interval %d (lambda_max < 0): [%.3f, %.3f]\n', r_pick, I_all(r_pick,1), I_all(r_pick,2));

% --- Lambda criterion: |lambda_max| >= n_decay / t_transient ---
lam_min = n_decay / t_transient;
i_int   = i_on(r_pick):i_off(r_pick);                  % grid points of the interval
ok_lam  = lle_row(i_int) <= -lam_min;
if ~any(ok_lam)
    error('No ks in the interval has |lambda_max| >= %.4f (n_decay = %g, t_transient = %g).', ...
        lam_min, n_decay, t_transient);
end
[j1, j2] = longestRun(ok_lam);
ks_a = ks_lle(i_int(j1));
ks_b = ks_lle(i_int(j2));
fprintf('Lambda criterion |lambda_max| >= %.4f: ks in [%.3f, %.3f]\n', lam_min, ks_a, ks_b);

% --- Pre-trim with the RMS file (use_rms_file) ---
rms_file = fullfile('mat', sprintf('RMS_vs_ks_%s_%s_T%d.mat', FMtype, topology_list{topology_idx}, T_end));   % this topology, same time profile
if use_rms_file
    if ~isfile(rms_file)
        error('%s not found (use_rms_file = true): run Motif_RMS_vs_ks_AllTopologies.m with time_profile = ''%s''.', rms_file, time_profile);
    end
    Rm = load(rms_file, 'RMS_topo', 'ks_vec', 'sync_tol', 'tspan', 't_transient');
    if isempty(rms_file_tol), rms_file_tol = Rm.sync_tol; end
    ks_rms  = Rm.ks_vec(:).';
    rms_row = Rm.RMS_topo(CasoBs, :);
    ok_rms  = ks_rms >= ks_a & ks_rms <= ks_b & rms_row < rms_file_tol;
    if ~any(ok_rms)
        error('No ks in [%.3f, %.3f] has rms(e) < %g in %s.', ks_a, ks_b, rms_file_tol, rms_file);
    end
    [j1, j2] = longestRun(ok_rms);
    fprintf('RMS file rms(e) < %g (its t_end = %g, transient = %g): ks in [%.3f, %.3f]\n', ...
        rms_file_tol, Rm.tspan(end), Rm.t_transient, ks_rms(j1), ks_rms(j2));
    ks_a = max(ks_a, ks_rms(j1));
    ks_b = min(ks_b, ks_rms(j2));
end

ks_sweep  = linspace(ks_a, ks_b, N_sweep);
lle_sweep = interp1(ks_lle, lle_row, ks_sweep);
fprintf('Sweep: ks in [%.3f, %.3f], %d values (lambda_max from %+.4f to %+.4f)\n', ...
    ks_a, ks_b, N_sweep, lle_sweep(1), lle_sweep(end));

% =========================================================================
% STEP B: INDICES, SCALE AND NODE-1 POOL
% =========================================================================
p_scale   = OscPar(1:3);
p_scale   = [p_scale p_scale];                      % a b c a b c
% Attractor scale x_ref of Eq. (eq_normalized_sensitivity), one value per
% oscillator (RMS of the nominal node trajectory, AttractorScale.m):
% sbar_ij = s_ij * q_j / x_ref
xref = AttractorScale(OsID);
err_idx   = 4:12;                                   % error rows of x
Ncoef     = numel(err_idx) * Nq;                    % 54 coefficients (q_j-major)
sens_rows = zeros(Ncoef, 1);
sens_scal = zeros(Ncoef, 1);
kk = 1;
for j = 1:Nq
    for i = err_idx
        sens_rows(kk) = Nx + (j-1)*Nx + i;
        sens_scal(kk) = p_scale(j) / xref;          % sbar_ij = s_ij * q_j / x_ref
        kk = kk + 1;
    end
end

xs_pool = sym('xs_pool', [n 1], 'real');
F_pool  = subs(F, [X; P.'], [xs_pool; OscPar(1:3).']);
f_pool  = matlabFunction(F_pool, 'Vars', {xs_pool});
rng(base_seed, 'twister');
[~, Xp] = ode45(@(t, x) f_pool(x), [0, pool_T_transient + (1:max_runs) * pool_spacing], ...
                rand(n, 1), odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
X1_pool = Xp(2:end, :);                             % max_runs x n

t_samp   = t_transient + (1:floor((T_end - t_transient) / dt_sample)) * dt_sample;
n_samp   = numel(t_samp);
% Blow-up test on node 1 AND the errors (rows 1:12): an escape of the
% trajectory stops the run at once instead of letting ode45 shrink its step
% until it fails; the solver warning of a failed run is not printed
warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
ode_opts = odeset('RelTol', ode_RelTol, 'AbsTol', ode_AbsTol, ...
                  'Events', @(t, Z) blowupEvent(Z, [1:n, err_idx], err_bound));

% Initial conditions (identical at every ks). Runs start ON the
% synchronization manifold: X1(0) from the pool (run r uses point r of the
% attractor), e(0) = 0 and S(0) = 0. With zero mismatches e stays exactly
% 0, so the sensitivity of node 1 (S_X1) does not enter s(e) (its coupling
% terms are proportional to e) and its rows are held at 0 (keep). In
% topologies where node 1 is driven, S_X1 grows like the chaos of node 1
% and, with e(0) ~= 0, the term S_X1*e made s(e) diverge (and ode45 crawl).
Z0_all = zeros(Nz, max_runs);
for r = 1:max_runs
    Z0_all(1:n, r) = X1_pool(r, :).';                % e(0) = 0, S(0) = 0
end
keep = ones(Nz, 1);
keep(setdiff(Nx+1:Nz, sens_rows)) = 0;              % rows of S of X1 held at 0
ic_mode = 'manifold';                               % stored in the .mat file

if replot_only
    % Redraw only: the saved profile replaces Steps C and D
    in_file = fullfile('mat', sprintf('SensitivityNorm_vs_ks_%s_%s_B%d_N%d_T%d.mat', ...
        FMtype, topology_list{topology_idx}, CasoBs, N_sweep, T_end));
    if ~isfile(in_file), error('replot_only: %s not found.', in_file); end
    if ~ismember('xref', who('-file', in_file))
        error('replot_only: %s has no xref field (run ShiftSensitivity_xref.m first).', in_file);
    end
    load(in_file, 'ks_sweep', 'lle_sweep', 'norm_dB', 'norm_unc', 'comp_dB', 'n_blown', 'no_influence', ...
        'S_worst', 'S_mean', 'ks_best', 'S_best', 'n_used', 'conv_ks', 'rms_e_med', 'ks_ok', ...
        'norm_dB_all', 'norm_unc_all', 'xref');
    ks_plot = [min(ks_sweep(ks_ok)), max(ks_sweep(ks_ok))];
    [~, k_best] = min(S_worst);
    fprintf('replot_only: profile loaded from %s\n', in_file);
else
% =========================================================================
% STEP C: RUNS AT EVERY ks
%   RMS_all(k, r, c): RMS of sbar of coefficient c, run r, ks index k
% =========================================================================
RMS_all = NaN(N_sweep, max_runs, Ncoef);
rms_e   = NaN(N_sweep, max_runs);                % rms(e) of each run in the RMS window
n_blown = zeros(1, N_sweep);
n_used  = zeros(1, N_sweep);                     % runs performed at each ks
conv_ks = false(1, N_sweep);                     % convergence test met at each ks
t_start = tic;
for k = 1:N_sweep
    Parval = [zeros(1, 6), ks_sweep(k)];            % ALL mismatches = 0
    f_ode  = @(t, Z) keep .* odeMotifSensitivityMismatch(t, Z, Parval);
    r_done   = 0;
    hist_r   = [];
    hist_idx = [];
    while r_done < max_runs && ~conv_ks(k)
        batch = r_done + 1 : min(r_done + batch_size, max_runs);
        for r = batch
            [t_out, Z_out, te] = ode45(f_ode, [0, t_samp], Z0_all(:, r), ode_opts);
            if ~isempty(te) || numel(t_out) < n_samp + 1 || any(~isfinite(Z_out(:)))
                n_blown(k) = n_blown(k) + 1;
                continue;
            end
            Sb = Z_out(2:end, sens_rows) .* sens_scal.';           % samples x 54
            RMS_all(k, r, :) = sqrt(mean(Sb.^2, 1));
            rms_e(k, r) = sqrt(mean(sum(Z_out(2:end, err_idx).^2, 2)));
        end
        r_done = batch(end);

        % Convergence test (same rule as SensitivityIndexMonteCarlo.m)
        Rk = reshape(RMS_all(k, 1:r_done, :), r_done, Ncoef);
        Rk = Rk(all(isfinite(Rk), 2), :);
        nv = size(Rk, 1);
        if nv < 2, continue; end
        act     = mean(Rk, 1) > zero_rel_tol * max(mean(Rk, 1));   % non-zero coefficients
        idx_now = NaN(1, Ncoef);
        unc_now = NaN(1, Ncoef);
        D = 20*log10(Rk(:, act));
        idx_now(act) = median(D, 1);
        unc_now(act) = 1.2533 * std(D, 0, 1) / sqrt(nv);
        hist_r   = [hist_r; r_done];                              %#ok<AGROW>
        hist_idx = [hist_idx; idx_now];                           %#ok<AGROW>
        i_ref = find(hist_r <= r_done - conv_window, 1, 'last');
        if nv >= min_runs && ~isempty(i_ref)
            delta   = max(abs(idx_now - hist_idx(i_ref, :)), [], 'omitnan');
            sem_max = max(unc_now, [], 'omitnan');
            conv_ks(k) = delta < conv_tol_dB && sem_max < conv_sem_dB;
        end
    end
    n_used(k) = r_done;
    fprintf('  ks = %.4f (%2d/%d): lambda_max = %+.4f, %3d runs (%s), %d blown up  (%.1f min)\n', ...
        ks_sweep(k), k, N_sweep, lle_sweep(k), r_done, ternaryText(conv_ks(k), 'converged', 'max_runs, NOT converged'), ...
        n_blown(k), toc(t_start)/60);
end

% =========================================================================
% STEP D: STATISTICS PER ks
%   comp_dB(k, c) : median over runs of 20 log10 RMS (-Inf = zero coefficient)
%   norm_dB(k, v) : 20 log10 || [s_x s_y s_z] || (Eq. 23 on linear medians)
%   norm_unc(k, v): 1.2533 * std(per-run norm in dB) / sqrt(n)
%   v = (j-1)*3 + p : mismatch q_j, error vector p (e12, e13, e23)
% =========================================================================
comp_dB  = NaN(N_sweep, Ncoef);
norm_dB  = NaN(N_sweep, Nq * num_pairs);
norm_unc = NaN(N_sweep, Nq * num_pairs);
for k = 1:N_sweep
    R  = reshape(RMS_all(k, :, :), max_runs, Ncoef); % max_runs x 54 (NaN beyond n_used)
    R  = R(all(isfinite(R), 2), :);                  % valid runs
    nv = size(R, 1);
    if nv < 2, continue; end
    is_zero = mean(R, 1) <= zero_rel_tol * max(mean(R, 1));
    cmp = median(20*log10(R), 1);
    cmp(is_zero) = -Inf;
    comp_dB(k, :) = cmp;
    for j = 1:Nq
        for p = 1:num_pairs
            v    = (j-1)*num_pairs + p;
            cols = (j-1)*numel(err_idx) + (p-1)*n + (1:n);
            if all(is_zero(cols)), continue; end    % no influence: stays NaN
            lin = 10.^(cmp(cols) / 20);
            norm_dB(k, v)  = 20*log10(norm(lin));
            D = 20*log10(sqrt(sum(R(:, cols).^2, 2)));
            norm_unc(k, v) = 1.2533 * std(D) / sqrt(nv);
        end
    end
end
% --- Internal double check: synchronization error of the same runs ---
rms_e_med = median(rms_e, 2, 'omitnan');          % N_sweep x 1
ks_ok     = true(N_sweep, 1);
if use_internal_rms_check
    ks_ok = rms_e_med <= rms_e_tol;
    if any(~ks_ok)
        fprintf('Flagged by the internal RMS check (median rms(e) > %g), left out of the figures:\n', rms_e_tol);
        fprintf('  ks = %.4f: median rms(e) = %.3g, lambda_max = %+.4f\n', ...
            [ks_sweep(~ks_ok); rms_e_med(~ks_ok).'; lle_sweep(~ks_ok)]);
    else
        fprintf('Internal RMS check: all %d ks values have median rms(e) <= %g.\n', N_sweep, rms_e_tol);
    end
    if ~any(ks_ok)
        error('Every ks was flagged by the internal RMS check (rms_e_tol = %g).', rms_e_tol);
    end
end
% Plotted values: flagged ks removed (the unfiltered ones are saved too)
norm_dB_all  = norm_dB;
norm_unc_all = norm_unc;
norm_dB(~ks_ok, :)  = NaN;
norm_unc(~ks_ok, :) = NaN;
ks_plot = [min(ks_sweep(ks_ok)), max(ks_sweep(ks_ok))];   % x-range of the valid points
no_influence = all(isnan(norm_dB), 1);             % vectors that are zero at every ks

% Most robust coupling: smallest worst-case sensitivity over the sweep
S_worst = max(norm_dB, [], 2, 'omitnan');
S_mean  = 20*log10(mean(10.^(norm_dB / 20), 2, 'omitnan'));
[S_best, k_best] = min(S_worst);
ks_best = ks_sweep(k_best);
fprintf('\nMost robust ks (smallest worst-case ||s||): ks = %.4f, worst case %.2f dB, lambda_max = %+.4f\n', ...
    ks_best, S_best, lle_sweep(k_best));
fprintf('Largest uncertainty of any norm: %.3f dB; runs excluded: %d of %d\n', ...
    max(norm_unc(:), [], 'omitnan'), sum(n_blown), sum(n_used));
fprintf('Runs per ks: %d to %d; not converged at %d of %d ks values\n', min(n_used), max(n_used), nnz(~conv_ks), N_sweep);

end   % if replot_only

% =========================================================================
% STEP E: FIGURE 1 - ||s_qj(e_ij)|| (dB) vs ks, one tile per mismatch
% =========================================================================
osc_names = {'R\"{o}ssler', 'Lorenz'};
par_tex   = {'\delta a_2', '\delta b_2', '\delta c_2', '\delta a_3', '\delta b_3', '\delta c_3'};
err_tex   = {'12', '13', '23'};
col_e     = [0.850 0.325 0.098; 0.494 0.184 0.556; 0.466 0.674 0.188];   % e12, e13, e23
mrk_e     = {'o', 's', '^'};
head_str  = sprintf('%s - %s: $\\mathbf{b}_%d \\to \\mathbf{c}_%d$', osc_names{OsID}, motif_names{topology_idx}, find(bs), find(cs));

fig1 = figure('Name', 'Sensitivity norm vs ks', 'Position', fig1_size);
tl1  = tiledlayout(fig1, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
y_all = norm_dB(isfinite(norm_dB));
y_lim = [floor(min(y_all) / 10) * 10, ceil(max(y_all) / 10) * 10];
ax1 = gobjects(1, Nq);
for j = 1:Nq
    ax = nexttile(tl1);
    ax1(j) = ax;
    hold(ax, 'on');
    h   = gobjects(1, num_pairs);
    lbl = cell(1, num_pairs);
    for p = 1:num_pairs
        v = (j-1)*num_pairs + p;
        lbl{p} = sprintf('$\\mathbf{e}_{%s}$', err_tex{p});
        if no_influence(v)
            continue;                               % Zero vector: no curve and no legend entry
        end
        if show_errorbars                           % error bars only, no markers
            errorbar(ax, ks_sweep, norm_dB(:, v), 2*norm_unc(:, v), 'LineStyle', 'none', ...
                'Color', col_e(p,:), 'LineWidth', 1.0, 'CapSize', 3, 'HandleVisibility', 'off');
        end
        % Markers every marker_step points, shifted by one point per curve
        % so the markers of e12, e13 and e23 never fall on the same ks
        mk_idx = p:marker_step:N_sweep;
        h(p) = plot(ax, ks_sweep, norm_dB(:, v), '-', 'Marker', mrk_e{p}, 'MarkerIndices', mk_idx, ...
            'Color', col_e(p,:), 'MarkerFaceColor', col_e(p,:), 'MarkerSize', 6, 'LineWidth', 1.5);
    end
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontSize', font_ticks, 'TickLabelInterpreter', 'latex');
    xlim(ax, [0, L.kf]);                            % Whole computed ks range [0, kf]
    ylim(ax, y_lim);
    title(ax, sprintf('$%s$', par_tex{j}), 'Interpreter', 'latex', 'FontSize', font_subtitle);
    if show_panel_letters
        text(ax, 0.97, 0.93, sprintf('(%s)', char('a' + j - 1)), 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'Interpreter', 'latex', 'FontSize', font_ticks);   % Panel letter
    end
    ok = isgraphics(h);
    legend(ax, h(ok), lbl(ok), 'Interpreter', 'latex', 'FontSize', font_legend, ...
        'Location', 'southoutside', 'Orientation', 'horizontal');   % Below the panel: never covers the curves
end
xlabel(tl1, '$k_s$', 'Interpreter', 'latex', 'FontSize', font_xlabel);
ylabel(tl1, 'RMS Sensitivity (dB)', 'Interpreter', 'latex', 'FontSize', font_ylabel);
if show_title
    title(tl1, sprintf('%s: $\\mathbf{b}_%d \\to \\mathbf{c}_%d$', motif_names{topology_idx}, find(bs), find(cs)), ...
        'Interpreter', 'latex', 'FontSize', font_title);
end

% =========================================================================
% STEP F: FIGURE 2 - ROBUSTNESS SUMMARY
%   (a) lambda_max(ks) over [0, kf]: swept interval and most robust ks
%   (b) worst-case and mean ||s|| vs ks
%   (c) ||s|| vs |lambda_max| (log-log): linear response ||s|| ~ 1/|lambda|
% =========================================================================
fig2 = figure('Name', 'Robustness summary');
tl2  = tiledlayout(fig2, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl2);
hold(ax, 'on');
y_lo = min(lle_row, [], 'omitnan');
y_hi = max(lle_row, [], 'omitnan');
patch(ax, [ks_a ks_b ks_b ks_a], [y_lo y_lo y_hi y_hi], [0.85 0.93 1.0], 'EdgeColor', 'none');
plot(ax, ks_lle, lle_row, 'b-', 'LineWidth', 1.5);
plot(ax, ks_sweep(ks_ok), lle_sweep(ks_ok), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
plot(ax, ks_sweep(~ks_ok), lle_sweep(~ks_ok), 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);   % flagged by the RMS check
yline(ax, 0, 'k--', 'LineWidth', 1.0);
yline(ax, -lam_min, 'k:', 'LineWidth', 1.0);     % lambda criterion
xline(ax, ks_best, 'r-', 'LineWidth', 1.5);
grid(ax, 'on');
set(ax, 'FontSize', font_ticks, 'TickLabelInterpreter', 'latex', 'Layer', 'top');
xlim(ax, [0, L.kf]);
ylim(ax, [y_lo, y_hi]);
xlabel(ax, '$k_s$', 'Interpreter', 'latex', 'FontSize', font_xlabel);
ylabel(ax, '$\lambda_{\max}$', 'Interpreter', 'latex', 'FontSize', font_ylabel);
title(ax, '(a) Swept interval', 'Interpreter', 'latex', 'FontSize', font_subtitle);

ax = nexttile(tl2);
hold(ax, 'on');
plot(ax, ks_sweep, S_worst, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'r');
plot(ax, ks_sweep, S_mean,  'k-s', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'k');
plot(ax, ks_best, S_best, 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'y');
grid(ax, 'on');
set(ax, 'FontSize', font_ticks, 'TickLabelInterpreter', 'latex');
xlim(ax, ks_plot + [-1 1] * 0.02 * max(diff(ks_plot), eps));
xlabel(ax, '$k_s$', 'Interpreter', 'latex', 'FontSize', font_xlabel);
ylabel(ax, 'Sensitivity (dB)', 'Interpreter', 'latex', 'FontSize', font_ylabel);
legend(ax, {'worst case $\max \|s_{q_j}(\mathbf{e}_{ij})\|$', 'mean over vectors', ...
    sprintf('most robust $k_s = %.2f$', ks_best)}, 'Interpreter', 'latex', 'FontSize', font_legend, 'Location', 'best');
title(ax, '(b) Robustness vs $k_s$', 'Interpreter', 'latex', 'FontSize', font_subtitle);

ax = nexttile(tl2);
hold(ax, 'on');
lam_dB = 20*log10(abs(lle_sweep(:)));
for j = 1:Nq
    for p = 1:num_pairs
        v = (j-1)*num_pairs + p;
        if no_influence(v), continue; end
        plot(ax, lam_dB, norm_dB(:, v), 'o', 'Color', col_e(p,:), 'MarkerSize', 4, 'MarkerFaceColor', col_e(p,:));
    end
end
x_ref = [min(lam_dB), max(lam_dB)];
y_ref = max(S_worst) - (x_ref - min(lam_dB));       % slope -1 in dB/dB: ||s|| ~ 1/|lambda|
h_ref = plot(ax, x_ref, y_ref, 'k--', 'LineWidth', 1.5);
grid(ax, 'on');
set(ax, 'FontSize', font_ticks, 'TickLabelInterpreter', 'latex');
xlabel(ax, '$20\log_{10}|\lambda_{\max}|$ (dB)', 'Interpreter', 'latex', 'FontSize', font_xlabel);
ylabel(ax, '$\|s_{q_j}(\mathbf{e}_{ij})\|$ (dB)', 'Interpreter', 'latex', 'FontSize', font_ylabel);
h_e = gobjects(1, num_pairs);
for p = 1:num_pairs
    h_e(p) = plot(ax, NaN, NaN, 'o', 'Color', col_e(p,:), 'MarkerFaceColor', col_e(p,:));
end
legend(ax, [h_e, h_ref], {'$\mathbf{e}_{12}$', '$\mathbf{e}_{13}$', '$\mathbf{e}_{23}$', ...
    'slope $-1$: $\|s\| \propto 1/|\lambda_{\max}|$'}, 'Interpreter', 'latex', 'FontSize', font_legend, 'Location', 'best');
title(ax, '(c) Sensitivity vs transverse stability', 'Interpreter', 'latex', 'FontSize', font_subtitle);
title(tl2, head_str, 'Interpreter', 'latex', 'FontSize', font_title);

% =========================================================================
% STEP G: SAVE (mat/ and eps/)
% =========================================================================
if ~exist('eps', 'dir'), mkdir('eps'); end
if ~exist('mat', 'dir'), mkdir('mat'); end
all_ax = findall([fig1, fig2], 'Type', 'axes');
for kk = 1:numel(all_ax)
    all_ax(kk).Toolbar.Visible = 'off';
end
file_tag = sprintf('%s_%s_B%d_N%d_T%d', FMtype, topology_list{topology_idx}, CasoBs, N_sweep, T_end);   % T in the name: profiles do not overwrite each other
exportgraphics(fig1, fullfile('eps', sprintf('SensitivityNorm_vs_ks_%s.eps', file_tag)), 'ContentType', 'vector');
exportgraphics(fig2, fullfile('eps', sprintf('SensitivityRobustness_vs_ks_%s.eps', file_tag)), 'ContentType', 'vector');
if ~replot_only
save(fullfile('mat', sprintf('SensitivityNorm_vs_ks_%s.mat', file_tag)), ...
    'ks_sweep', 'lle_sweep', 'norm_dB', 'norm_unc', 'comp_dB', 'RMS_all', 'n_blown', ...
    'no_influence', 'S_worst', 'S_mean', 'ks_best', 'S_best', 'I_all', 'r_pick', 'ks_a', 'ks_b', ...
    'OsID', 'FMtype', 'OscPar', 'topology_idx', 'CasoBs', 'max_runs', 'n_used', 'conv_ks', 'batch_size', 'min_runs', ...
    'conv_window', 'conv_tol_dB', 'conv_sem_dB', 'T_end', 't_transient', ...
    'dt_sample', 'base_seed', 'X1_pool', 'p_scale', 'time_profile', 'n_decay', 'lam_min', 'use_rms_file', 'rms_file_tol', 'use_internal_rms_check', 'rms_e_tol', 'rms_e', 'rms_e_med', 'ks_ok', ...
    'norm_dB_all', 'norm_unc_all', 'lle_file', 'ic_mode', 'xref');
end
if replot_only
    fprintf('Saved: eps/SensitivityNorm_vs_ks_%s.eps, eps/SensitivityRobustness_vs_ks_%s.eps (mat file not rewritten)\n', ...
        file_tag, file_tag);
else
    fprintf('Saved: eps/SensitivityNorm_vs_ks_%s.eps, eps/SensitivityRobustness_vs_ks_%s.eps, mat/SensitivityNorm_vs_ks_%s.mat\n', ...
        file_tag, file_tag, file_tag);
end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function [value, isterminal, direction] = blowupEvent(Z, err_idx, err_bound)
    % Stops a run when any |e_ij| exceeds err_bound or the state is NaN/Inf
    if any(~isfinite(Z))
        value = -1;
    else
        value = err_bound - max(abs(Z(err_idx)));
    end
    isterminal = 1;
    direction  = 0;
end

function [i1, i2] = longestRun(mask)
    % First and last index of the longest run of true values in mask
    d  = diff([0, mask(:).', 0]);
    on = find(d == 1);
    of = find(d == -1) - 1;
    [~, r] = max(of - on);
    i1 = on(r);
    i2 = of(r);
end

function str = ternaryText(cond, str_true, str_false)
    if cond
        str = str_true;
    else
        str = str_false;
    end
end
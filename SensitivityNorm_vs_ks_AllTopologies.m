% =========================================================================
% SCRIPT: SensitivityNorm_vs_ks_AllTopologies.m
% DESCRIPTION: Batch version of SensitivityNorm_vs_ks.m. Computes the RMS
%              sensitivity vector index ||s_qj(e_ij)|| (dB) against ks, for
%              every topology in topology_run and every scheme in
%              scheme_run that has a VALID synchronization interval in
%              mat/LLE_vs_ks_<Osc>_AllTopologies.mat (lambda_max < 0 over at
%              least min_width, also after the lambda criterion); the RMS
%              file only trims its edges. The valid cases and their swept
%              intervals are printed before any run. The method is the
%              same as SensitivityNorm_vs_ks.m (same settings, trimming,
%              convergence test, seeds and statistics).
%              No figures. Each (topology, scheme) is saved in its own
%              file as soon as it finishes,
%                  mat/SensitivityNorm_vs_ks_<Osc>_<topology>_B<C>_N<N>_T<T_end>.mat
%              (same name and variables as SensitivityNorm_vs_ks.m, so a
%              file written by that script is reused). A stopped run keeps
%              every finished combination and, with skip_existing = true,
%              continues where it stopped. A scheme without a valid
%              interval is not computed and gets NO file (a file of such a
%              scheme left by an earlier run is deleted when
%              delete_invalid_files = true).
%              Per combination the sensitivity equations are generated in
%              odeSensNorm_<topology>_B<C>.m (deleted afterwards unless
%              keep_ode_files = true).
%              Runs start on the synchronization manifold (e(0) = 0,
%              S(0) = 0, X1(0) from a pool of points of the attractor); the
%              sensitivity of node 1 is held at 0 (it does not enter s(e)
%              when e = 0). Files of the earlier method (random e(0), no
%              ic_mode field) are recomputed and overwritten.
%              With use_parallel = true the 25 ks values of a combination
%              run in parallel (parfor over ks); runs and seeds are the same
%              as in serial, so the results do not depend on the mode.
%              At the end SensitivityNorm_Consolidate.m joins all the files
%              (consolidate_at_end).
% =========================================================================
clear; clc; close all;

% =========================================================================
% CONFIGURATION
% =========================================================================
OsID = 2;                       % 1=Rossler, 2=Lorenz

% Topologies and schemes computed in this run
%   1 chain, 2 auxiliar, 3 relay, 4 mixed_master, 5 competitive,
%   6 mixed_slave, 7 ring_bidir_chain, 8 ring_bidir_double_mixed,
%   9 ring_bidir_mixed, 10 ring_unidir_chain, 11 ring_bidir_auxiliar,
%   12 ring_unidir_auxiliar, 13 ring_bidir_competitive
topology_run  = 8:13;
scheme_run    = 1:9;            % CasoBs (b_i -> c_j, SISObc)
skip_existing = true;           % true: a combination whose file exists is not recomputed
% Files of the earlier method (random e(0), no ic_mode field): false = keep
% them and skip (nothing already computed is overwritten); true = recompute
% and overwrite them with the current method (e(0) = 0)
recompute_old_method = false;
% Only schemes with a valid synchronization interval in the LLE file get a
% .mat file. A file of an invalid scheme (left by an earlier run) is deleted
% when delete_invalid_files = true.
delete_invalid_files = true;

% --- ks sweep inside a synchronization interval of the LLE file ---
N_sweep       = 25;             % Number of ks values
interval_pick_table = [0 1];    % per OsID: 0 = longest interval; k = k-th interval (in ks order)
                                % Lorenz: 1 = first window (b3 -> c3 has a second window at large ks)
interval_pick = interval_pick_table(OsID);
ks_max_table  = [Inf 20];       % per OsID: LLE intervals cut at ks_max (Lorenz table: ks in [0, 20])
% Excluded combinations [topology_idx CasoBs] per OsID: no first window in [0, ks_max]
exclude_table = {zeros(0, 2), [4 9]};   % Lorenz: linear mixed master, b3 -> c3
% Edge trimming (lambda criterion): |lambda_max| * t_transient >= n_decay
n_decay       = 5;              % 5 -> |lambda_max| >= 5 / t_transient
% Valid scheme (decided from the LLE file): it has an interval with
% lambda_max < 0 at least min_width wide, and after the lambda criterion it
% is still at least min_width wide. Narrower pieces (a few grid points where
% lambda_max dips below 0, e.g. periodic windows or finite-time noise) are
% not synchronization. Invalid schemes are not computed and get no file.
min_width     = 0.1;            % minimum width of a valid interval (ks units)
% Internal double check: a ks is flagged when the median over runs of
% rms(e) for t >= t_transient exceeds rms_e_tol, e = [e12; e13; e23]
use_internal_rms_check = true;
rms_e_tol              = 1e-3;
% Pre-trim with mat/RMS_vs_ks_<Osc>_<topology>_T<T_end>.mat (Motif_RMS_vs_ks_AllTopologies.m
% run with the SAME time_profile). A topology whose RMS file does not exist
% yet is skipped (listed at the start) and can be computed in a later run.
use_rms_file = true;
rms_file_tol = [];              % [] = sync_tol stored in the RMS file

% --- Runs per ks: same convergence test as SensitivityIndexMonteCarlo.m ---
% Safety cap chosen from the oscillator (only reached if a ks never converges):
%              Rossler  Lorenz  Chen   HR    Duffing  vanderPol
max_runs_table = [ 100,  200,    200,   200,  200,     200 ];
max_runs     = max_runs_table(OsID);
batch_size   = 10;              % Runs between convergence checks
min_runs     = 10;              % Valid runs before the test is applied
conv_window  = 10;              % Compare with the index conv_window runs earlier
conv_tol_dB  = 1;               % max |change of the index| over the window (dB)
conv_sem_dB  = 1;               % max uncertainty of the index (dB)
% Time horizon profile (sets t_transient and T_end):
%   'short'  : t_transient =  200, T_end = 1000  -> lambda_min = 0.025
%   'medium' : t_transient =  500, T_end = 2500  -> lambda_min = 0.010
%   'long'   : t_transient = 1000, T_end = 5000  -> lambda_min = 0.005
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
base_seed   = 42;               % Run r uses rng(base_seed + r) (same at every ks and combination)
pool_T_transient = 500;         % Node-1 pool: discarded transient
pool_spacing     = 50;          % Node-1 pool: time between stored points
err_bound    = 1e3;             % Blow-up bound on |X1| and |e_ij| (run excluded)
zero_rel_tol = 1e-6;            % Zero coefficient: mean RMS <= zero_rel_tol * max at that ks

% --- Computing ---
use_parallel   = true;          % true: parfor over ks (Parallel Computing Toolbox), tested on the first combination
n_workers      = [];            % [] = default of the cluster profile; or a number of workers
keep_ode_files = false;         % false: delete odeSensNorm_<topology>_B<C>.m after each combination
show_progress  = true;          % true: a progress line after every batch of runs of every ks (live, also in parallel)

% At the end, join all the files with SensitivityNorm_Consolidate.m (it
% only writes the unified file when all 13 x 9 files exist, unless its
% allow_incomplete = true)
consolidate_at_end = true;

% =========================================================================
% STEP 1: OSCILLATOR, TOPOLOGIES, SCHEMES AND LLE FILE
% =========================================================================
num_oscillators = 3;
syms a b c ks t real
[F, OscPar, FMtype] = OscillatorType(OsID);
syms x(t) y(t) z(t)
X = [x(t); y(t); z(t)];
P = [a b c];
n = length(X);                  % State dimension per node (3)
num_pairs = 3;                  % e12, e13, e23
Nx = n * (1 + num_pairs);       % 12 states [X1; e12; e13; e23]
Nq = 6;                         % 6 mismatches [da2 db2 dc2 da3 db3 dc3]
Nz = Nx * (1 + Nq);             % 84 states [x; vec(S)]

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
edges_list = { ...
    [1 2; 2 3], ...                                 % chain
    [1 2; 1 3], ...                                 % auxiliar
    [1 2; 2 1; 2 3; 3 2], ...                       % relay
    [1 2; 2 3; 3 2], ...                            % mixed_master
    [1 2; 3 2], ...                                 % competitive
    [2 1; 2 3; 3 2], ...                            % mixed_slave
    [1 2; 2 1; 2 3; 3 2; 3 1; 1 3], ...             % ring_bidir_chain
    [1 2; 2 3; 3 2; 3 1; 1 3], ...                  % ring_bidir_double_mixed
    [1 2; 2 3; 3 2; 3 1], ...                       % ring_bidir_mixed
    [1 2; 2 3; 3 1], ...                            % ring_unidir_chain
    [1 2; 1 3; 2 3; 3 2], ...                       % ring_bidir_auxiliar
    [1 2; 1 3; 2 3], ...                            % ring_unidir_auxiliar
    [1 2; 2 1; 2 3; 1 3]};                          % ring_bidir_competitive
num_topo = numel(topology_list);
if any(topology_run < 1 | topology_run > num_topo) || any(scheme_run < 1 | scheme_run > 9)
    error('topology_run must be in 1..%d and scheme_run in 1..9.', num_topo);
end

scheme_bc    = zeros(9, 2);     % [i j] of b_i -> c_j
scheme_names = strings(1, 9);
for CasoBs = 1:9
    [bs, cs] = SISObc(CasoBs, n);
    scheme_bc(CasoBs, :)  = [find(bs), find(cs)];
    scheme_names(CasoBs)  = sprintf('b%d_to_c%d', find(bs), find(cs));
end

lle_file = fullfile('mat', sprintf('LLE_vs_ks_%s_AllTopologies.mat', FMtype));
if ~isfile(lle_file)
    error('%s not found (run MotifLLE_vs_ks_AllTopologies.m / LLE_Consolidate.m first).', lle_file);
end
L = load(lle_file, 'LLE_array', 'ks_vec', 'kf');
ks_lle  = L.ks_vec(:).';
% Cut at ks_max and drop the excluded combinations (NaN = not synchronized)
L.LLE_array(:, :, ks_lle > ks_max_table(OsID)) = NaN;
ex = exclude_table{OsID};
for r = 1:size(ex, 1)
    L.LLE_array(ex(r,1), ex(r,2), :) = NaN;
end

% Topologies whose RMS file does not exist yet are skipped in this run
if use_rms_file
    has_rms = false(size(topology_run));
    for it = 1:numel(topology_run)
        rms_file = fullfile('mat', sprintf('RMS_vs_ks_%s_%s_T%d.mat', FMtype, topology_list{topology_run(it)}, T_end));
        has_rms(it) = isfile(rms_file);
    end
    if ~any(has_rms)
        error('No RMS file found for the topologies of this run (run Motif_RMS_vs_ks_AllTopologies.m with time_profile = ''%s'').', ...
            time_profile);
    end
    if ~all(has_rms)
        fprintf(['RMS file not found yet (time_profile ''%s''); these topologies are skipped in this run:\n  %s\n', ...
            'Computing: %s\n'], time_profile, strjoin(topology_list(topology_run(~has_rms)), ', '), ...
            strjoin(topology_list(topology_run(has_rms)), ', '));
        topology_run = topology_run(has_rms);
    end
end
lam_min = n_decay / t_transient;
if ~exist('mat', 'dir'), mkdir('mat'); end

% =========================================================================
% STEP 2: INDICES, SCALE, NODE-1 POOL AND INITIAL CONDITIONS
%   Identical for every combination (common random numbers)
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

t_samp = t_transient + (1:floor((T_end - t_transient) / dt_sample)) * dt_sample;

% Runs start ON the synchronization manifold: X1(0) from the pool (run r
% uses point r of the attractor), e(0) = 0 and S(0) = 0. With zero
% mismatches e stays exactly 0, so the sensitivity of node 1 (S_X1) does
% not enter s(e) (its coupling terms are proportional to e) and its rows
% are held at 0 in runOneKs. In topologies where node 1 is driven, S_X1
% grows like the chaos of node 1 and, with e(0) ~= 0, the term S_X1*e made
% s(e) diverge (and ode45 crawl).
Z0_all = zeros(Nz, max_runs);
for r = 1:max_runs
    Z0_all(1:n, r) = X1_pool(r, :).';                % e(0) = 0, S(0) = 0
end

% Settings passed to the run and statistics functions
cfg = struct('t_samp', t_samp, 'n_samp', numel(t_samp), 'RelTol', ode_RelTol, 'AbsTol', ode_AbsTol, ...
    'err_idx', err_idx, 'err_bound', err_bound, 'sens_rows', sens_rows, 'sens_scal', sens_scal, ...
    'Ncoef', Ncoef, 'Nq', Nq, 'n', n, 'num_pairs', num_pairs, 'max_runs', max_runs, ...
    'batch_size', batch_size, 'min_runs', min_runs, 'conv_window', conv_window, ...
    'conv_tol_dB', conv_tol_dB, 'conv_sem_dB', conv_sem_dB, 'zero_rel_tol', zero_rel_tol, ...
    'dt_sample', dt_sample);

% =========================================================================
% STEP 3: PARALLEL POOL (the serial = parfor test runs on the first
%         combination, once its ODE file exists)
% =========================================================================
pool_tested = false;
if use_parallel
    fprintf('Starting the parallel pool...\n');
    par_ok = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
    if ~par_ok
        warning('Parallel Computing Toolbox not available: running in serial.');
    else
        try
            pool = gcp('nocreate');
            if isempty(pool)
                if isempty(n_workers)
                    pool = parpool;
                else
                    pool = parpool(n_workers);
                end
            elseif ~isempty(n_workers) && pool.NumWorkers ~= n_workers
                delete(pool);
                pool = parpool(n_workers);
            end
            fprintf('  Pool running with %d workers.\n', pool.NumWorkers);
        catch ME
            warning('Parallel pool could not be started (%s): running in serial.', ME.message);
            par_ok = false;
        end
    end
    use_parallel = par_ok;
end

% =========================================================================
% STEP 3b: VALID SCHEMES AND INTERVALS FROM THE LLE FILE (before any run)
%   LLE    : interval with lambda_max < 0, trimmed by the lambda criterion
%   swept  : the LLE interval trimmed at the edges by the RMS file
%   Only the cases with a swept interval are computed and saved.
% =========================================================================
fprintf('\n%s, profile %s (t_transient = %g, T_end = %g), N_sweep = %d, lambda_min = %.4f, min_width = %g\n', ...
    FMtype, time_profile, t_transient, T_end, N_sweep, lam_min, min_width);
fprintf('\nSwept interval of each case (-- = no valid interval in the LLE file, RMS = removed by the RMS file):\n');
fprintf('%-32s', 'Topology');
fprintf(' %15s', string(compose('b%d->c%d', scheme_bc(scheme_run, :))));
fprintf('\n');
n_valid = 0;
for topology_idx = topology_run
    lle_topo_p = reshape(L.LLE_array(topology_idx, :, :), 9, []);
    if use_rms_file
        Rp = load(fullfile('mat', sprintf('RMS_vs_ks_%s_%s_T%d.mat', FMtype, topology_list{topology_idx}, T_end)), ...
            'RMS_topo', 'ks_vec', 'sync_tol');
        rms_tol_p = rms_file_tol;
        if isempty(rms_tol_p), rms_tol_p = Rp.sync_tol; end
    end
    fprintf('%-32s', motif_names{topology_idx});
    for CasoBs = scheme_run
        LI = lleInterval(ks_lle, lle_topo_p(CasoBs, :), interval_pick, lam_min, min_width, false);
        if ~isempty(LI.reason)
            fprintf(' %15s', '--');
            continue;
        end
        ks_a_p = LI.ks_a;
        ks_b_p = LI.ks_b;
        if use_rms_file
            [ks_a_p, ks_b_p, reason_p] = rmsTrim(ks_a_p, ks_b_p, Rp.ks_vec(:).', Rp.RMS_topo(CasoBs, :), ...
                rms_tol_p, min_width, false);
            if ~isempty(reason_p)
                fprintf(' %15s', 'RMS');
                continue;
            end
        end
        fprintf(' %15s', sprintf('[%.2f, %.2f]', ks_a_p, ks_b_p));
        n_valid = n_valid + 1;
    end
    fprintf('\n');
end
fprintf('%d of %d cases have a valid interval and will be computed (unless their file exists).\n', ...
    n_valid, numel(topology_run) * numel(scheme_run));

% =========================================================================
% STEP 4: SWEEP OVER TOPOLOGIES AND SCHEMES
% =========================================================================
t_total = tic;
n_computed = 0;
for topology_idx = topology_run
    topo_name = topology_list{topology_idx};
    edges     = edges_list{topology_idx};
    fprintf('\n=== Topology %d: %s ===\n', topology_idx, motif_names{topology_idx});
    lle_topo = reshape(L.LLE_array(topology_idx, :, :), 9, []);
    if all(isnan(lle_topo(:)))
        fprintf('  No LLE data for this topology in %s: skipped.\n', lle_file);
        continue;
    end

    % RMS file of this topology (same time profile)
    rms_file = fullfile('mat', sprintf('RMS_vs_ks_%s_%s_T%d.mat', FMtype, topo_name, T_end));
    rms_tol_used = rms_file_tol;
    if use_rms_file
        Rm = load(rms_file, 'RMS_topo', 'ks_vec', 'sync_tol', 'tspan', 't_transient');
        if isempty(rms_tol_used), rms_tol_used = Rm.sync_tol; end
    end

    for CasoBs = scheme_run
        bc       = scheme_bc(CasoBs, :);
        file_tag = sprintf('%s_%s_B%d_N%d_T%d', FMtype, topo_name, CasoBs, N_sweep, T_end);
        out_file = fullfile('mat', sprintf('SensitivityNorm_vs_ks_%s.mat', file_tag));
        fprintf('\n  %s, b_%d -> c_%d (CasoBs = %d)\n', motif_names{topology_idx}, bc(1), bc(2), CasoBs);
        t_combo = tic;

        % --- STEP A1: valid scheme and interval from the LLE file ---
        lle_row = lle_topo(CasoBs, :);
        LI = lleInterval(ks_lle, lle_row, interval_pick, lam_min, min_width, true);
        if ~isempty(LI.reason)
            fprintf('    Not computed: %s\n', LI.reason);
            if isfile(out_file) && delete_invalid_files
                delete(out_file);
                fprintf('    Deleted %s (scheme without a valid interval).\n', out_file);
            end
            continue;
        end

        % --- Existing file ---
        if isfile(out_file)
            Sx = struct();
            vars_file = who('-file', out_file);
            if ismember('ic_mode', vars_file) && ismember('norm_dB', vars_file) && ~ismember('xref', vars_file)
                error(['%s was computed without the attractor scale x_ref ' ...
                    '(no xref field): run ShiftSensitivity_xref.m first.'], out_file);
            end
            if ismember('status', vars_file)        % files of SensitivityNorm_vs_ks.m have no status
                Sx = load(out_file, 'status');
            end
            if isfield(Sx, 'status') && strcmp(Sx.status, 'no_interval')
                delete(out_file);                   % marker of an earlier version: recompute
                fprintf('    Deleted the old no_interval file %s: computing.\n', out_file);
            elseif ~ismember('ic_mode', vars_file) && recompute_old_method   % random e(0): recompute
                fprintf('    %s is of the earlier method (random e(0)): recomputing (overwritten when saved).\n', out_file);
            elseif skip_existing
                if ~ismember('ic_mode', vars_file)
                    fprintf('    %s exists (earlier method, random e(0); recompute_old_method = false): skipped.\n', out_file);
                    continue;
                end
                fprintf('    %s exists: skipped.\n', out_file);
                continue;
            end
        end

        % --- STEP A2: edge trimming with the RMS file ---
        SW = struct('I_all', LI.I_all, 'r_pick', LI.r_pick, 'ks_a', LI.ks_a, 'ks_b', LI.ks_b);
        if use_rms_file
            [SW.ks_a, SW.ks_b, reason] = rmsTrim(LI.ks_a, LI.ks_b, Rm.ks_vec(:).', Rm.RMS_topo(CasoBs, :), ...
                rms_tol_used, min_width, true);
            if ~isempty(reason)
                warning('%s, b_%d -> c_%d: valid in the LLE file but %s: not computed.', ...
                    motif_names{topology_idx}, bc(1), bc(2), reason);
                continue;
            end
        end

        % Common fields of the output file
        S_out = struct('OsID', OsID, 'FMtype', FMtype, 'OscPar', OscPar, ...
            'topology_idx', topology_idx, 'topology_name', topo_name, ...
            'motif_name', motif_names{topology_idx}, 'edges', edges, ...
            'CasoBs', CasoBs, 'scheme_bc', bc, 'scheme_name', char(scheme_names(CasoBs)), ...
            'N_sweep', N_sweep, 'interval_pick', interval_pick, 'max_runs', max_runs, ...
            'batch_size', batch_size, 'min_runs', min_runs, 'conv_window', conv_window, ...
            'conv_tol_dB', conv_tol_dB, 'conv_sem_dB', conv_sem_dB, 'T_end', T_end, ...
            't_transient', t_transient, 'dt_sample', dt_sample, 'base_seed', base_seed, ...
            'X1_pool', X1_pool, 'p_scale', p_scale, 'xref', xref, 'time_profile', time_profile, ...
            'n_decay', n_decay, 'lam_min', lam_min, 'use_rms_file', use_rms_file, ...
            'rms_file_tol', rms_tol_used, 'use_internal_rms_check', use_internal_rms_check, ...
            'rms_e_tol', rms_e_tol, 'err_bound', err_bound, 'zero_rel_tol', zero_rel_tol, ...
            'ode_RelTol', ode_RelTol, 'ode_AbsTol', ode_AbsTol, 'lle_file', lle_file, ...
            'min_width', min_width);
        S_out.ic_mode  = 'manifold';                % e(0) = 0, S(0) = 0, S_X1 held at 0
        S_out.rms_file = '';
        if use_rms_file, S_out.rms_file = rms_file; end
        S_out.I_all  = SW.I_all;
        S_out.r_pick = SW.r_pick;
        S_out.ks_a   = SW.ks_a;
        S_out.ks_b   = SW.ks_b;

        ks_sweep  = linspace(SW.ks_a, SW.ks_b, N_sweep);
        lle_sweep = interp1(ks_lle, lle_row, ks_sweep);
        fprintf('    Sweep: ks in [%.3f, %.3f], %d values (lambda_max from %+.4f to %+.4f)\n', ...
            SW.ks_a, SW.ks_b, N_sweep, lle_sweep(1), lle_sweep(end));

        % --- Sensitivity equations of this topology and scheme ---
        ode_name = sprintf('odeSensNorm_%s_B%d', topo_name, CasoBs);
        fprintf('    Generating %s.m ...\n', ode_name);
        buildSensitivityODE(OsID, edges, CasoBs, ode_name);
        clear(ode_name);
        rehash;
        if use_parallel
            refreshWorkers(ode_name);               % Workers pick up the new file
        end
        if use_parallel && ~pool_tested
            use_parallel = poolTest(ode_name, ks_sweep(ceil(N_sweep/2)), Z0_all(:, 1:2), cfg);
            pool_tested  = true;
        end
        fprintf('    Computing mode: %s\n', ternaryText(use_parallel, 'parallel (parfor over ks)', 'serial'));

        % --- STEP C: runs at every ks ---
        out = cell(1, N_sweep);
        % Live progress: in parallel the workers send their lines to the
        % client through a DataQueue; in serial they are printed directly
        q = [];
        if use_parallel && show_progress
            q = parallel.pool.DataQueue;
            afterEach(q, @(msg) fprintf('%s', msg));
        end
        fprintf('    Running %d ks values (a line after every %d runs of each ks)...\n', N_sweep, batch_size);
        if use_parallel
            parfor k = 1:N_sweep
                out{k} = runOneKs(ode_name, ks_sweep(k), Z0_all, cfg, ...
                    struct('k', k, 'N', N_sweep, 'q', q, 'on', show_progress));
            end
        else
            for k = 1:N_sweep
                out{k} = runOneKs(ode_name, ks_sweep(k), Z0_all, cfg, ...
                    struct('k', k, 'N', N_sweep, 'q', q, 'on', show_progress));
            end
        end
        fprintf('    Summary of the %d ks values:\n', N_sweep);
        RMS_all = NaN(N_sweep, max_runs, Ncoef);
        rms_e   = NaN(N_sweep, max_runs);
        n_blown = zeros(1, N_sweep);
        n_used  = zeros(1, N_sweep);
        conv_ks = false(1, N_sweep);
        for k = 1:N_sweep
            RMS_all(k, :, :) = reshape(out{k}.RMS, 1, max_runs, Ncoef);
            rms_e(k, :) = out{k}.rms_e;
            n_blown(k)  = out{k}.n_blown;
            n_used(k)   = out{k}.n_used;
            conv_ks(k)  = out{k}.conv;
            fprintf('    ks = %.4f (%2d/%d): lambda_max = %+.4f, %3d runs (%s), %d blown up  (%.1f min)\n', ...
                ks_sweep(k), k, N_sweep, lle_sweep(k), n_used(k), ...
                ternaryText(conv_ks(k), 'converged', 'max_runs, NOT converged'), n_blown(k), out{k}.time_min);
        end

        % --- STEP D: statistics per ks (same as SensitivityNorm_vs_ks.m) ---
        [comp_dB, norm_dB, norm_unc] = sweepStats(RMS_all, cfg);
        rms_e_med = median(rms_e, 2, 'omitnan');
        ks_ok     = true(N_sweep, 1);
        if use_internal_rms_check
            ks_ok = rms_e_med <= rms_e_tol;
            if any(~ks_ok)
                fprintf('    Flagged by the internal RMS check (median rms(e) > %g):\n', rms_e_tol);
                fprintf('      ks = %.4f: median rms(e) = %.3g, lambda_max = %+.4f\n', ...
                    [ks_sweep(~ks_ok); rms_e_med(~ks_ok).'; lle_sweep(~ks_ok)]);
            end
        end
        norm_dB_all  = norm_dB;
        norm_unc_all = norm_unc;
        norm_dB(~ks_ok, :)  = NaN;
        norm_unc(~ks_ok, :) = NaN;
        no_influence = all(isnan(norm_dB), 1);
        S_worst = max(norm_dB, [], 2, 'omitnan');
        S_mean  = 20*log10(mean(10.^(norm_dB / 20), 2, 'omitnan'));
        if any(ks_ok)
            status  = 'ok';
            reason  = '';
            ks_plot = [min(ks_sweep(ks_ok)), max(ks_sweep(ks_ok))];
            [S_best, k_best] = min(S_worst);
            ks_best = ks_sweep(k_best);
        else
            status  = 'all_flagged';
            reason  = sprintf('every ks flagged by the internal RMS check (rms_e_tol = %g)', rms_e_tol);
            ks_plot = [NaN NaN];
            S_best  = NaN;
            k_best  = NaN;
            ks_best = NaN;
            fprintf('    %s\n', reason);
        end

        % --- STEP E: save this combination ---
        S_out.status       = status;
        S_out.reason       = reason;
        S_out.ks_sweep     = ks_sweep;
        S_out.lle_sweep    = lle_sweep;
        S_out.norm_dB      = norm_dB;
        S_out.norm_unc     = norm_unc;
        S_out.norm_dB_all  = norm_dB_all;
        S_out.norm_unc_all = norm_unc_all;
        S_out.comp_dB      = comp_dB;
        S_out.RMS_all      = RMS_all;
        S_out.rms_e        = rms_e;
        S_out.rms_e_med    = rms_e_med;
        S_out.ks_ok        = ks_ok;
        S_out.no_influence = no_influence;
        S_out.S_worst      = S_worst;
        S_out.S_mean       = S_mean;
        S_out.S_best       = S_best;
        S_out.k_best       = k_best;
        S_out.ks_best      = ks_best;
        S_out.ks_plot      = ks_plot;
        S_out.n_blown      = n_blown;
        S_out.n_used       = n_used;
        S_out.conv_ks      = conv_ks;
        S_out.run_time_min  = toc(t_combo) / 60;
        S_out.date_computed = datestr(now, 'yyyy-mm-dd HH:MM');
        save(out_file, '-struct', 'S_out');
        n_computed = n_computed + 1;
        fprintf('    Most robust ks = %.4f (worst case %.2f dB); not converged at %d of %d ks; %d flagged; %.1f min\n', ...
            ks_best, S_best, nnz(~conv_ks), N_sweep, nnz(~ks_ok), S_out.run_time_min);
        fprintf('    Saved: %s\n', out_file);

        if ~keep_ode_files
            clear(ode_name);
            if isfile([ode_name '.m']), delete([ode_name '.m']); end
        end
    end
end
fprintf('\n%d combinations computed in %.1f h\n', n_computed, toc(t_total)/3600);

% =========================================================================
% STEP 5: UNIFIED FILE (SensitivityNorm_Consolidate.m, same settings)
% =========================================================================
if consolidate_at_end
    if exist('SensitivityNorm_Consolidate.m', 'file')
        fprintf('\nConsolidating the per-combination files...\n');
        consolidate_OsID    = OsID;
        consolidate_rule    = struct('interval_pick', interval_pick, 'n_decay', n_decay, ...
            'min_width', min_width, 'use_rms_file', use_rms_file, 'rms_file_tol', rms_file_tol);
        consolidate_N_sweep = N_sweep;
        consolidate_T_end   = T_end;
        run('SensitivityNorm_Consolidate.m');
    else
        warning('SensitivityNorm_Consolidate.m not found in the path: unified file not written.');
    end
end

% =========================================================================
% LOCAL FUNCTION: VALID INTERVAL FROM THE LLE FILE
%   Intervals with lambda_max < 0 at least min_width wide; the longest (or
%   the interval_pick-th) is trimmed by the lambda criterion
%   |lambda_max| >= lam_min and must still be at least min_width wide.
%   LI.reason = '' when the scheme is valid.
% =========================================================================
function LI = lleInterval(ks_lle, lle_row, interval_pick, lam_min, min_width, verbose)
    LI = struct('I_all', zeros(0, 2), 'r_pick', NaN, 'ks_a', NaN, 'ks_b', NaN, 'reason', '');
    if all(isnan(lle_row))
        LI.reason = 'no LLE data';
        return;
    end
    mask  = lle_row < 0;                            % NaN (blow-up) counts as not synchronized
    d     = diff([0, mask, 0]);
    i_on  = find(d == 1);
    i_off = find(d == -1) - 1;
    keep  = ks_lle(i_off) - ks_lle(i_on) >= min_width;   % drop pieces narrower than min_width
    i_on  = i_on(keep);
    i_off = i_off(keep);
    if isempty(i_on)
        LI.reason = sprintf('no interval with lambda_max < 0 at least %g wide', min_width);
        return;
    end
    LI.I_all = [ks_lle(i_on).', ks_lle(i_off).'];
    if verbose, fprintf('    LLE intervals (lambda_max < 0):%s\n', sprintf(' [%.3f, %.3f]', LI.I_all.')); end
    if interval_pick == 0
        [~, r_pick] = max(LI.I_all(:,2) - LI.I_all(:,1));
    elseif interval_pick <= size(LI.I_all, 1)
        r_pick = interval_pick;
    else
        LI.reason = sprintf('interval_pick = %d but only %d intervals', interval_pick, size(LI.I_all, 1));
        return;
    end
    LI.r_pick = r_pick;
    i_int  = i_on(r_pick):i_off(r_pick);
    ok_lam = lle_row(i_int) <= -lam_min;
    if ~any(ok_lam)
        LI.reason = sprintf('no ks in [%.3f, %.3f] has |lambda_max| >= %.4f', LI.I_all(r_pick,1), LI.I_all(r_pick,2), lam_min);
        return;
    end
    [j1, j2] = longestRun(ok_lam);
    LI.ks_a = ks_lle(i_int(j1));
    LI.ks_b = ks_lle(i_int(j2));
    if LI.ks_b - LI.ks_a < min_width
        LI.reason = sprintf('after the lambda criterion [%.3f, %.3f] is narrower than %g', LI.ks_a, LI.ks_b, min_width);
        return;
    end
    if verbose, fprintf('    Interval %d; lambda criterion |lambda_max| >= %.4f: ks in [%.3f, %.3f]\n', r_pick, lam_min, LI.ks_a, LI.ks_b); end
end

% =========================================================================
% LOCAL FUNCTION: EDGE TRIMMING WITH THE RMS FILE
%   Keeps the longest run of ks in [ks_a, ks_b] with rms(e) < rms_tol; the
%   result must be at least min_width wide (reason = '' when it is).
% =========================================================================
function [ks_a, ks_b, reason] = rmsTrim(ks_a, ks_b, ks_rms, rms_row, rms_tol, min_width, verbose)
    reason = '';
    ok_rms = ks_rms >= ks_a & ks_rms <= ks_b & rms_row < rms_tol;
    if ~any(ok_rms)
        reason = sprintf('no ks in [%.3f, %.3f] has rms(e) < %g in the RMS file', ks_a, ks_b, rms_tol);
        return;
    end
    [j1, j2] = longestRun(ok_rms);
    if verbose, fprintf('    RMS file rms(e) < %g: ks in [%.3f, %.3f]\n', rms_tol, ks_rms(j1), ks_rms(j2)); end
    ks_a = max(ks_a, ks_rms(j1));
    ks_b = min(ks_b, ks_rms(j2));
    if ks_b - ks_a < min_width
        reason = sprintf('the RMS file leaves [%.3f, %.3f], narrower than %g', ks_a, ks_b, min_width);
    end
end

% =========================================================================
% LOCAL FUNCTION: RUNS AT ONE ks UNTIL THE CONVERGENCE TEST IS MET
%   (same loop as STEP C of SensitivityNorm_vs_ks.m)
% =========================================================================
function R = runOneKs(ode_name, ks_val, Z0_all, cfg, prog)
    % prog: k, N (ks index and count), q (DataQueue or []), on (print progress)
    t_ks     = tic;
    f_sens   = str2func(ode_name);
    Parval   = [zeros(1, 6), ks_val];               % ALL mismatches = 0
    % On the synchronization manifold (e = 0) the sensitivity of node 1 does
    % not enter s(e): its rows are held at 0 (exact), so its chaotic growth
    % cannot slow ode45
    keep = ones(size(Z0_all, 1), 1);
    keep(setdiff(cfg.n + numel(cfg.err_idx) + 1 : size(Z0_all, 1), cfg.sens_rows)) = 0;
    f_ode    = @(t, Z) keep .* f_sens(t, Z, Parval);
    % Blow-up test on node 1 AND the errors (rows 1:12): an escape of the
    % trajectory from the attractor stops the run at once instead of letting
    % ode45 shrink its step until it fails (hours per run). Such runs are
    % excluded as before; the solver warning of a failed run is not printed.
    warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
    ode_opts = odeset('RelTol', cfg.RelTol, 'AbsTol', cfg.AbsTol, ...
                      'Events', @(t, Z) blowupEvent(Z, [1:cfg.n, cfg.err_idx], cfg.err_bound));
    RMS      = NaN(cfg.max_runs, cfg.Ncoef);
    rms_e    = NaN(1, cfg.max_runs);
    n_blown  = 0;
    conv     = false;
    r_done   = 0;
    hist_r   = [];
    hist_idx = [];
    while r_done < cfg.max_runs && ~conv
        delta   = NaN;
        sem_max = NaN;
        batch = r_done + 1 : min(r_done + cfg.batch_size, cfg.max_runs);
        for r = batch
            [t_out, Z_out, te] = ode45(f_ode, [0, cfg.t_samp], Z0_all(:, r), ode_opts);
            if ~isempty(te) || numel(t_out) < cfg.n_samp + 1 || any(~isfinite(Z_out(:)))
                n_blown = n_blown + 1;
                continue;
            end
            Sb = Z_out(2:end, cfg.sens_rows) .* cfg.sens_scal.';           % samples x 54
            RMS(r, :) = sqrt(mean(Sb.^2, 1));
            rms_e(r)  = sqrt(mean(sum(Z_out(2:end, cfg.err_idx).^2, 2)));
        end
        r_done = batch(end);

        % Convergence test (same rule as SensitivityIndexMonteCarlo.m)
        Rk = RMS(1:r_done, :);
        Rk = Rk(all(isfinite(Rk), 2), :);
        nv = size(Rk, 1);
        if nv < 2, continue; end
        act     = mean(Rk, 1) > cfg.zero_rel_tol * max(mean(Rk, 1));   % non-zero coefficients
        idx_now = NaN(1, cfg.Ncoef);
        unc_now = NaN(1, cfg.Ncoef);
        D = 20*log10(Rk(:, act));
        idx_now(act) = median(D, 1);
        unc_now(act) = 1.2533 * std(D, 0, 1) / sqrt(nv);
        hist_r   = [hist_r; r_done];                                  %#ok<AGROW>
        hist_idx = [hist_idx; idx_now];                               %#ok<AGROW>
        i_ref = find(hist_r <= r_done - cfg.conv_window, 1, 'last');
        if nv >= cfg.min_runs && ~isempty(i_ref)
            delta   = max(abs(idx_now - hist_idx(i_ref, :)), [], 'omitnan');
            sem_max = max(unc_now, [], 'omitnan');
            conv    = delta < cfg.conv_tol_dB && sem_max < cfg.conv_sem_dB;
        end
        report(prog, sprintf(['      ks %2d/%d = %.4f: %3d runs, %d blown up, change %5.2f dB, ', ...
            'uncertainty %5.2f dB (%.1f min)\n'], prog.k, prog.N, ks_val, r_done, n_blown, ...
            delta, sem_max, toc(t_ks) / 60));
    end
    report(prog, sprintf('    > ks %2d/%d = %.4f finished: %d runs, %s, %d blown up (%.1f min)\n', ...
        prog.k, prog.N, ks_val, r_done, ternaryText(conv, 'converged', 'NOT converged'), n_blown, toc(t_ks) / 60));
    R = struct('RMS', RMS, 'rms_e', rms_e, 'n_blown', n_blown, 'n_used', r_done, ...
               'conv', conv, 'time_min', toc(t_ks) / 60);
end

function report(prog, msg)
    % Progress line: printed in serial, sent to the client in parallel
    if ~prog.on, return; end
    if isempty(prog.q)
        fprintf('%s', msg);
    else
        send(prog.q, msg);
    end
end

% =========================================================================
% LOCAL FUNCTION: STATISTICS PER ks (same as STEP D of SensitivityNorm_vs_ks.m)
%   comp_dB(k, c) : median over runs of 20 log10 RMS (-Inf = zero coefficient)
%   norm_dB(k, v) : 20 log10 || [s_x s_y s_z] || (Eq. 23 on linear medians)
%   norm_unc(k, v): 1.2533 * std(per-run norm in dB) / sqrt(n)
%   v = (j-1)*3 + p : mismatch q_j, error vector p (e12, e13, e23)
% =========================================================================
function [comp_dB, norm_dB, norm_unc] = sweepStats(RMS_all, cfg)
    N_sweep  = size(RMS_all, 1);
    n_err    = numel(cfg.err_idx);
    comp_dB  = NaN(N_sweep, cfg.Ncoef);
    norm_dB  = NaN(N_sweep, cfg.Nq * cfg.num_pairs);
    norm_unc = NaN(N_sweep, cfg.Nq * cfg.num_pairs);
    for k = 1:N_sweep
        R  = reshape(RMS_all(k, :, :), cfg.max_runs, cfg.Ncoef);
        R  = R(all(isfinite(R), 2), :);             % valid runs
        nv = size(R, 1);
        if nv < 2, continue; end
        is_zero = mean(R, 1) <= cfg.zero_rel_tol * max(mean(R, 1));
        cmp = median(20*log10(R), 1);
        cmp(is_zero) = -Inf;
        comp_dB(k, :) = cmp;
        for j = 1:cfg.Nq
            for p = 1:cfg.num_pairs
                v    = (j-1)*cfg.num_pairs + p;
                cols = (j-1)*n_err + (p-1)*cfg.n + (1:cfg.n);
                if all(is_zero(cols)), continue; end    % no influence: stays NaN
                lin = 10.^(cmp(cols) / 20);
                norm_dB(k, v)  = 20*log10(norm(lin));
                D = 20*log10(sqrt(sum(R(:, cols).^2, 2)));
                norm_unc(k, v) = 1.2533 * std(D) / sqrt(nv);
            end
        end
    end
end

% =========================================================================
% LOCAL FUNCTION: SENSITIVITY EQUATIONS -> <ode_name>.m
%   Same construction as PART 1 of SensitivityNorm_vs_ks.m:
%   Z = [x (12); vec(S) (72)], x = [X1; e12; e13; e23],
%   Parval = [da2 db2 dc2 da3 db3 dc3 ks], a, b, c fixed at OscPar.
% =========================================================================
function buildSensitivityODE(OsID, edges, CasoBs, ode_name)
    num_oscillators = 3;
    syms a b c ks t real
    [F, OscPar, ~] = OscillatorType(OsID);
    syms x(t) y(t) z(t)
    X = [x(t); y(t); z(t)];
    n = length(X);

    % Parameter mismatches of nodes 2 and 3 (node 1 nominal)
    syms da2 db2 dc2
    syms da3 db3 dc3
    P   = [a b c];
    DP2 = [da2 db2 dc2];
    DP3 = [da3 db3 dc3];
    DP_all = {sym(zeros(1, numel(P))), DP2, DP3};
    Par = [P DP2 DP3 ks];
    num_pairs = (num_oscillators * (num_oscillators - 1)) / 2;
    num_edges = size(edges, 1);
    [bs, cs] = SISObc(CasoBs, n);

    % States and vector fields of the three nodes
    XAll = [];
    FX_All = [];
    for k = 1:num_oscillators
        eval(sprintf('syms x%d(t) y%d(t) z%d(t)', k, k, k));
        eval(sprintf('X%d = [x%d(t); y%d(t); z%d(t)];', k, k, k, k));
        eval(sprintf('FX%d = subs(F, X, X%d);', k, k));
        eval(sprintf('FX%d = subs(FX%d, P, P + DP_all{%d});', k, k, k));
        XAll   = [XAll;   eval(sprintf('X%d', k))];   %#ok<AGROW>
        FX_All = [FX_All; eval(sprintf('FX%d', k))];  %#ok<AGROW>
    end
    Xj  = XAll;
    FXj = FX_All;

    idx = cell(num_oscillators, 1);
    for k = 1:num_oscillators
        idx{k} = (k-1)*n + (1:n);
    end

    % Error symbols
    E_array  = cell(num_pairs, 1);
    pair_ij  = zeros(num_pairs, 2);
    vars_E   = [];
    pair_idx = 1;
    for i = 1:num_oscillators
        for j = (i+1):num_oscillators
            pair_str = sprintf('e%d%d', i, j);
            pair_ij(pair_idx, :) = [i, j];
            eval(sprintf('syms %sx(t) %sy(t) %sz(t)', pair_str, pair_str, pair_str));
            eval(sprintf('%s = [%sx(t); %sy(t); %sz(t)];', pair_str, pair_str, pair_str, pair_str));
            eval(sprintf('E_array{%d} = %s;', pair_idx, pair_str));
            eval(sprintf('vars_E = [vars_E; %s];', pair_str));
            pair_idx = pair_idx + 1;
        end
    end

    % Coupling inputs in error variables
    Pidx = zeros(num_oscillators);
    for p = 1:num_pairs
        Pidx(pair_ij(p,1), pair_ij(p,2)) = p;
    end
    errDiff = @(s, k) (s < k) * E_array{Pidx(min(s,k), max(s,k))} ...
                    - (s > k) * E_array{Pidx(min(s,k), max(s,k))};   % X_s - X_k
    U = cell(num_oscillators, 1);
    for k = 1:num_oscillators
        U{k} = sym(zeros(n, 1));
    end
    for e = 1:num_edges
        src = edges(e,1);
        tgt = edges(e,2);
        U{tgt} = U{tgt} + ks * (bs * cs) * errDiff(src, tgt);
    end

    % Error equations: de_ij/dt = F_i(X_i) - F_j(X_i - e_ij) + u_i - u_j
    F_matrix_E = [];
    for p = 1:num_pairs
        i = pair_ij(p, 1);
        j = pair_ij(p, 2);
        Xi = Xj(idx{i});
        Fj = subs(F, P, P + DP_all{j});
        FXi_minus = subs(Fj, X, Xi - E_array{p});
        dEij = (FXj(idx{i}) - FXi_minus) + U{i} - U{j};
        F_matrix_E = [F_matrix_E; dEij];            %#ok<AGROW>
    end
    F_matrix_E = expand(F_matrix_E);
    for r = 1:length(F_matrix_E)
        F_matrix_E(r) = collect(F_matrix_E(r), Par);
    end
    F_matrix_E(4:6) = expand(subs(F_matrix_E(4:6), E_array{3}, E_array{2} - E_array{1}));
    for r = 4:6
        F_matrix_E(r) = collect(F_matrix_E(r), Par);
    end

    % Augmented system [X1; e12; e13; e23], X2 = X1 - e12
    X1v   = Xj(idx{1});
    F_X1  = FXj(idx{1}) + U{1};
    F_aug = [F_X1; F_matrix_E];
    F_aug = subs(F_aug, Xj(idx{2}), X1v - E_array{1});

    % Independent system [X1; e12; e13], nominal a, b, c
    F_red    = subs(F_aug(1:3*n), E_array{3}, E_array{2} - E_array{1});
    F_red    = subs(F_red, P, OscPar(1:3));
    vars_red = [X1v; E_array{1}; E_array{2}];
    if any(has(F_red, [E_array{3}; Xj(idx{2}); Xj(idx{3})]))
        error('Reduced system still depends on e23, X2 or X3.');
    end
    Par_run = [DP2 DP3 ks];

    % Sensitivity equations dS/dt = F_x S + F_q (e23 part without feedback)
    q   = [DP2 DP3];
    NxI = numel(vars_red);
    Nx  = NxI + n;
    Nq  = numel(q);
    syms s(t) [Nx Nq]
    S   = s(t);
    S_I = S(1:NxI, :);
    F_x   = jacobian(F_red, vars_red);
    F_q   = jacobian(F_red, q);
    F_S_I = F_x * S_I + F_q;
    i12 = n + (1:n);
    i13 = 2*n + (1:n);
    F_e23  = F_red(i13) - F_red(i12);
    F_S_23 = F_S_I(i13, :) - F_S_I(i12, :);
    F_x_all = [F_red; F_e23];
    F_S     = [F_S_I; F_S_23];
    vars_x  = [vars_red; E_array{3}];

    vars_sens = [vars_x; S(:)];
    F_sens    = [F_x_all; F_S(:)];
    odeFunction(F_sens, vars_sens, Par_run, 'File', [ode_name '.m'], 'Optimize', true);
end

% =========================================================================
% LOCAL FUNCTION: PARALLEL TEST ON THE FIRST COMBINATION
%   Two runs on a short horizon (t <= 20), serial vs parfor, must agree.
%   Returns false (serial computing) on any failure.
% =========================================================================
function ok = poolTest(ode_name, ks_test, Z0, cfg)
    fprintf('    Testing the parallel pool (serial vs parfor, 2 short runs)...\n');
    ok = true;
    try
        t_test = 0:cfg.dt_sample:20;
        Z_ser  = cell(1, 2);
        Z_par  = cell(1, 2);
        for r = 1:2
            Z_ser{r} = testRun(ode_name, ks_test, Z0(:, r), t_test, cfg);
        end
        t_par = tic;
        parfor r = 1:2
            Z_par{r} = testRun(ode_name, ks_test, Z0(:, r), t_test, cfg);
        end
        d = max(cellfun(@(A, B) max(abs(A(:) - B(:))), Z_ser, Z_par));
        if ~(d <= 1e-8)
            w_info = fetchOutputs(parfevalOnAll(@workerInfo, 1, ode_name), 'UniformOutput', false);
            fprintf('    client : %s\n', which(ode_name));
            fprintf('    workers: %s\n', strjoin(unique(w_info), ' || '));
            warning('parfor and serial results differ (%.2e): running in serial.', d);
            ok = false;
        else
            fprintf('    parfor test passed (serial = parfor, %.1f s).\n', toc(t_par));
        end
    catch ME
        warning('Parallel pool test failed (%s): running in serial.', ME.message);
        ok = false;
    end
end

function Z = testRun(ode_name, ks_val, Z0, t_test, cfg)
    f_sens = str2func(ode_name);
    Parval = [zeros(1, 6), ks_val];
    [~, Z] = ode45(@(t, Z) f_sens(t, Z, Parval), t_test, Z0, ...
                   odeset('RelTol', cfg.RelTol, 'AbsTol', cfg.AbsTol));
end

% =========================================================================
% LOCAL FUNCTION: WORKERS USE THE CLIENT FOLDER AND THE NEW ODE FILE
% =========================================================================
function refreshWorkers(ode_name)
    if isempty(gcp('nocreate'))
        return;
    end
    here = strrep(pwd, '''', '''''');
    pctRunOnAll(['cd(''' here ''')']);
    pctRunOnAll(['clear ' ode_name]);
    pctRunOnAll('rehash');
end

function info = workerInfo(ode_name)
    r    = rng;
    info = sprintf('%s (folder %s, rng %s)', which(ode_name), pwd, r.Type);
end

% =========================================================================
% LOCAL FUNCTIONS: BLOW-UP EVENT, LONGEST RUN, TEXT CHOICE
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

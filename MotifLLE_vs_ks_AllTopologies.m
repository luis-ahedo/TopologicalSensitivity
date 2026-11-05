% =========================================================================
% SCRIPT: MotifLLE_vs_ks_AllTopologies.m
% DESCRIPTION: Largest transverse Lyapunov exponent LLE(ks) for ks in
%              I_ks = [0, kf], for every coupling scheme CasoBs = 1 to 9
%              and the topologies listed in topology_run (all 13 by
%              default). No figures: each topology is saved in its own
%              .mat file as soon as it finishes,
%                  mat/LLE_vs_ks_<Osc>_<topology>.mat,
%              so the 13 topologies can be computed in groups (e.g. one
%              group per night) and a stopped run keeps every finished
%              topology. LLE_PlotAndTable.m reads all these files.
%              Per scheme the variational equations on e = 0 are rebuilt
%              and odeMotifVariational.m is regenerated; the exponents
%              come from MaxTransLyaExp.
%              Edge refinement (refine_edges): near the zero crossings of
%              LLE(ks) the slope of LLE(ks) is small, so the statistical
%              error of the finite-time average (about +/-0.005 for
%              Rossler, tf = 2000) moves the crossing by up to ~0.02 in
%              ks. The ks points around each crossing are recomputed with
%              a longer averaging time tf_refine; the window is widened
%              until the refined crossing lies inside it. refine_only =
%              true refines existing per-topology files without repeating
%              the coarse sweep.
% =========================================================================
clear; clc; close all;

% =========================================================================
% STEP 1: CONFIGURATION
% =========================================================================
OsID = 1;                       % 1=Rossler, 2=Lorenz
num_oscillators = 3;

% ks grid: I_ks = [0, kf] with step dks, kf chosen from the oscillator
%              Rossler  Lorenz  Chen    HR     Duffing  vanderPol
kf_table   = [ 5,       25,     25,     25,    25,      25    ];
dks        = 0.01;                          % ks step (grid resolution)
kf         = kf_table(OsID);
N_ks       = round(kf / dks) + 1;           % 501 (Rossler), 2501 (Lorenz)
ks_vec     = linspace(0, kf, N_ks);

% RK4 step and averaging time chosen from the oscillator. Checked for the
% chain b1 -> c1 against dt = 0.001: the step error is below the
% statistical error of the finite-time average, which is about +/-0.005
% (Rossler, tf = 2000) and +/-0.01 (Lorenz, tf = 1000), enough to place the
% zero crossings of LLE(ks) well within the ks grid spacing.
%              Rossler  Lorenz  Chen    HR     Duffing  vanderPol
tf_table   = [ 2000,    1000,   2000,   2000,  2000,    2000  ];
dt_table   = [ 0.01,    0.005,  0.005,  0.005, 0.005,   0.005 ];

opts = struct();
opts.tf          = tf_table(OsID);   % Maximum averaging time
opts.dt          = dt_table(OsID);   % RK4 fixed step
opts.T_transient = 100;         % Transient of node 1
opts.num_LE      = 1;           % Largest exponent only
opts.seed        = 42;          % Same perturbation for every ks, scheme, topology
opts.early_exit  = false;       % false: always integrate the full tf (true: stop early, see MaxTransLyaExp)

% --- Refinement of the synchronization edges (zero crossings of LLE) ---
%              Rossler  Lorenz  Chen    HR     Duffing  vanderPol
tf_refine_table = [ 20000,   10000,  20000,  20000, 20000,   20000 ];
refine_edges      = true;       % true: recompute the ks points around each zero crossing with tf_refine
refine_only       = false;      % true: load the existing per-topology files and ONLY refine their edges
                                %       (no coarse sweep); topologies already refined with the same
                                %       tf_refine are skipped
refine_half       = 0.05;       % Half-width (ks units) of the refined window around each crossing
refine_max_expand = 10;         % Max extra half-widths when the refined crossing leaves the window
opts_refine    = opts;
opts_refine.tf = tf_refine_table(OsID);

use_parallel = true;            % true: parfor over ks (Parallel Computing Toolbox), tested in STEP 3b
n_workers    = [];              % [] = default of the cluster profile; or a number of workers

% Topologies computed in this run (indices 1..13 of topology_list below),
% e.g. 1:4 tonight, 5:9 tomorrow, 10:13 the next night.
%   1 chain, 2 auxiliar, 3 relay, 4 mixed_master, 5 competitive,
%   6 mixed_slave, 7 ring_bidir_chain, 8 ring_bidir_double_mixed,
%   9 ring_bidir_mixed, 10 ring_unidir_chain, 11 ring_bidir_auxiliar,
%   12 ring_unidir_auxiliar, 13 ring_bidir_competitive
topology_run = 1:13;

% At the end of the run, join the per-topology files into
% mat/LLE_vs_ks_<Osc>_AllTopologies.mat with LLE_Consolidate.m (it only
% writes the file when all 13 topologies exist, unless its
% allow_incomplete = true)
consolidate_at_end = true;

% =========================================================================
% STEP 2: OSCILLATOR AND TOPOLOGY
% =========================================================================
syms a b c ks t real
[F, OscPar, FMtype] = OscillatorType(OsID);
syms x(t) y(t) z(t)
X   = [x(t); y(t); z(t)];
Par = [a, b, c, ks];
n   = length(X);
total_states = n * num_oscillators;
dim_Y = n * (num_oscillators - 1);    % 6 independent perturbations [de12; de13]

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
num_topo = numel(topology_list);

osc_names = {'R\"{o}ssler', 'Lorenz'};
motif_names = {'Linear chain', 'Linear auxiliary', 'Linear relay', ...
    'Linear mixed master', 'Linear competitive', 'Linear mixed slave', ...
    'Ring bidirectional chain', 'Ring bidirectional double mixed', ...
    'Ring bidirectional mixed', 'Ring unidirectional chain', ...
    'Ring bidirectional auxiliary', 'Ring unidirectional auxiliary', ...
    'Ring bidirectional competitive'};

% Node 1 initial condition (same seed as SimulateMotif.m)
rng(42);
initial_conditions = rand(total_states, 1);
opts.x0 = initial_conditions(1:n);

% =========================================================================
% STEP 3: OUTPUT FILES (one per topology)
%   mat/LLE_vs_ks_<Osc>_<topology>.mat contains
%   LLE_topo(CasoBs, ks index) = LLE(ks), 9 x N_ks
%   LLE_data.<bi_to_cj>        = LLE(ks), 1 x N_ks
% =========================================================================
if any(topology_run < 1 | topology_run > num_topo)
    error('topology_run must contain indices between 1 and %d.', num_topo);
end
scheme_names = strings(1, 9);
for CasoBs = 1:9
    [bs, cs] = SISObc(CasoBs, n);
    scheme_names(CasoBs) = sprintf('b%d_to_c%d', find(bs), find(cs));
end

mat_folder = 'mat';
if ~exist(mat_folder, 'dir')
    mkdir(mat_folder);
end

% =========================================================================
% STEP 3b: PARALLEL POOL TEST (only when use_parallel = true)
%   1) Parallel Computing Toolbox installed and licensed;
%   2) a pool is running or can be started (n_workers);
%   3) the workers see the generated odeMotifVariational.m and give the
%      same exponents as a serial run: two ks values of the first
%      topology and scheme with a short run (T_transient = 5, tf = 10).
%   Any failure prints the reason and falls back to serial computing.
% =========================================================================
if use_parallel
    fprintf('\nTesting the parallel pool...\n');
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

            % Same results serial vs parfor on a short run
            [bs, cs] = SISObc(1, n);
            buildVariationalODE(F, X, n, eval(['edges_' topology_list{topology_run(1)}]), bs, cs, Par);
            clear odeMotifVariational;
            rehash;
            refreshWorkers();                           % Workers: same folder, fresh file
            % Short run (chaos would amplify any round-off difference
            % between client and workers on long runs)
            opts_test             = opts;
            opts_test.T_transient = 5;
            opts_test.tf          = 10;
            opts_test.early_exit  = false;
            ks_test  = [ks_vec(2), ks_vec(end)];
            lle_ser  = zeros(1, 2);
            lle_par  = zeros(1, 2);
            for k_idx = 1:2
                lle_ser(k_idx) = MaxTransLyaExp([OscPar, ks_test(k_idx)], dim_Y, n, opts_test);
            end
            t_par = tic;
            parfor k_idx = 1:2
                lle_par(k_idx) = MaxTransLyaExp([OscPar, ks_test(k_idx)], dim_Y, n, opts_test);
            end
            if max(abs(lle_par - lle_ser)) > 1e-8
                % Diagnostics: which ODE file and folder each side uses
                fprintf('  serial LLE = [%s], parfor LLE = [%s]\n', num2str(lle_ser, '%.5f '), num2str(lle_par, '%.5f '));
                fprintf('  client : %s\n', which('odeMotifVariational'));
                w_info = fetchOutputs(parfevalOnAll(@workerInfo, 1), 'UniformOutput', false);
                fprintf('  workers: %s\n', strjoin(unique(w_info), ' || '));
                warning('parfor and serial results differ (%.2e): running in serial.', ...
                    max(abs(lle_par - lle_ser)));
                par_ok = false;
            else
                fprintf('  parfor test passed (serial = parfor, %.1f s).\n', toc(t_par));
            end
        catch ME
            warning('Parallel pool test failed (%s): running in serial.', ME.message);
            par_ok = false;
        end
    end
    use_parallel = par_ok;
end
fprintf('Computing mode: %s\n', ternaryText(use_parallel, 'parallel (parfor over ks)', 'serial'));

% =========================================================================
% STEP 4: SWEEP OVER TOPOLOGIES, SCHEMES AND ks
% =========================================================================
t_total = tic;
for i_run = 1:numel(topology_run)
    topology_idx = topology_run(i_run);
    edges_input = topology_list{topology_idx};
    edges       = eval(['edges_' edges_input]);
    fprintf('\n=== Topology %d (%d of %d in this run): %s ===\n', topology_idx, i_run, ...
        numel(topology_run), motif_names{topology_idx});
    t_topo = tic;
    LLE_data = struct();
    LLE_topo = NaN(9, N_ks);
    LLE_coarse   = NaN(9, N_ks);                    % Coarse sweep (tf = opts.tf)
    refined_mask = false(9, N_ks);                  % ks points recomputed with tf_refine
    filename_mat = fullfile(mat_folder, sprintf('LLE_vs_ks_%s_%s.mat', FMtype, edges_input));
    if refine_only
        if ~isfile(filename_mat)
            warning('%s not found: topology skipped (refine_only = true).', filename_mat);
            continue;
        end
        T_old = load(filename_mat);
        if isfield(T_old, 'opts_refine') && isfield(T_old, 'refined_mask') && ...
                T_old.opts_refine.tf == opts_refine.tf && any(T_old.refined_mask(:))
            fprintf('  Already refined with tf = %g: skipped.\n', opts_refine.tf);
            continue;
        end
        if isfield(T_old, 'LLE_coarse')
            LLE_coarse = T_old.LLE_coarse;
        else
            LLE_coarse = T_old.LLE_topo;
        end
        if ~isequal(size(LLE_coarse), [9 N_ks])
            warning('%s has another ks grid: topology skipped.', filename_mat);
            continue;
        end
    end

    for CasoBs = 1:9
        [bs, cs] = SISObc(CasoBs, n);
        fprintf('  CasoBs = %d / 9 (b_%d -> c_%d)\n', CasoBs, find(bs), find(cs));
        if refine_only && ~any(diff(LLE_coarse(CasoBs, :) < 0) ~= 0)
            LLE_data.(char(scheme_names(CasoBs))) = LLE_coarse(CasoBs, :);   % No crossing: nothing to refine
            LLE_topo(CasoBs, :) = LLE_coarse(CasoBs, :);
            continue;
        end

        % Variational equations on e = 0 -> odeMotifVariational.m
        buildVariationalODE(F, X, n, edges, bs, cs, Par);
        clear odeMotifVariational;
        rehash;
        if use_parallel
            refreshWorkers();                       % Workers pick up the new file
        end

        if refine_only
            LLE_row = LLE_coarse(CasoBs, :);           % Coarse sweep from the existing file
        else
            LLE_row = NaN(1, N_ks);
            if use_parallel
                parfor k_idx = 1:N_ks
                    LLE_row(k_idx) = MaxTransLyaExp([OscPar, ks_vec(k_idx)], dim_Y, n, opts);
                end
            else
                for k_idx = 1:N_ks
                    LLE_row(k_idx) = MaxTransLyaExp([OscPar, ks_vec(k_idx)], dim_Y, n, opts);
                end
            end
            LLE_coarse(CasoBs, :) = LLE_row;
        end

        % Edge refinement: longer average around each zero crossing
        if refine_edges
            [LLE_row, ref_mask] = refineEdges(LLE_row, ks_vec, OscPar, dim_Y, n, opts_refine, ...
                refine_half, refine_max_expand, use_parallel);
            refined_mask(CasoBs, :) = ref_mask;
            if any(ref_mask)
                fprintf('    refined %d ks points (tf = %g)\n', nnz(ref_mask), opts_refine.tf);
            end
        end

        LLE_data.(char(scheme_names(CasoBs))) = LLE_row;
        LLE_topo(CasoBs, :) = LLE_row;
    end

    run_time_min = toc(t_topo) / 60;
    fprintf('  Topology finished in %.1f min\n', run_time_min);

    % One file per topology, written as soon as it finishes
    topology_name = edges_input;
    motif_name    = motif_names{topology_idx};
    date_computed = datestr(now, 'yyyy-mm-dd HH:MM');
    save(filename_mat, 'LLE_topo', 'LLE_data', 'ks_vec', 'kf', 'N_ks', ...
        'scheme_names', 'topology_idx', 'topology_name', 'motif_name', 'edges', ...
        'osc_names', 'FMtype', 'OsID', 'OscPar', 'opts', 'run_time_min', 'date_computed', ...
        'LLE_coarse', 'refined_mask', 'refine_edges', 'opts_refine', 'refine_half');
    fprintf('  Saved: %s\n', filename_mat);
end
fprintf('\nTopologies %s finished in %.1f min\n', mat2str(topology_run), toc(t_total)/60);

% =========================================================================
% STEP 5: CONSOLIDATED FILE (LLE_Consolidate.m, same oscillator)
% =========================================================================
if consolidate_at_end
    if exist('LLE_Consolidate.m', 'file')
        fprintf('\nConsolidating the per-topology files...\n');
        consolidate_OsID = OsID;
        run('LLE_Consolidate.m');
    else
        warning('LLE_Consolidate.m not found in the path: consolidated file not written.');
    end
end

% =========================================================================
% LOCAL FUNCTION: VARIATIONAL EQUATIONS ON e = 0 -> odeMotifVariational.m
%   Error system (same formulation as MotifErrorSystem.m):
%     de_ij/dt = F(X_i) - F(X_i - e_ij) + u_i - u_j,  X_2 = X_1 - e12,
%     u_k = ks * sum_{s -> k} B*C*(X_s - X_k), each difference = +/- e,
%     e23 = e13 - e12 substituted in de12/dt and de13/dt.
%   Variational system (independent errors only, dY = [de12; de13]):
%     dx1/dt = F(x1),  d(dY)/dt = J_E(x1)|_{e=0} * dY.
% =========================================================================
function buildVariationalODE(F, X, n, edges, bs, cs, Par)
    pairs = [1 2; 1 3; 2 3];
    num_pairs = size(pairs, 1);
    dim_Y = n * num_pairs;

    xs = sym('xs', [n 1], 'real');                  % node 1 states
    es = sym('es', [dim_Y 1], 'real');              % errors [e12; e13; e23]
    dim_V = n * (num_pairs - 1);                    % independent errors [e12; e13]
    dY = sym('dY', [dim_V 1], 'real');              % perturbations [de12; de13]
    E  = cell(num_pairs, 1);
    for p = 1:num_pairs
        E{p} = es((p-1)*n + (1:n));
    end
    P = zeros(3);
    for p = 1:num_pairs
        P(pairs(p,1), pairs(p,2)) = p;
    end
    errDiff = @(s, k) (s < k) * E{P(min(s,k), max(s,k))} ...
                    - (s > k) * E{P(min(s,k), max(s,k))};   % X_s - X_k
    Fx = @(v) subs(F, X, v);                        % F evaluated at v

    % Coupling inputs in error variables
    U = {sym(zeros(n,1)); sym(zeros(n,1)); sym(zeros(n,1))};
    for e = 1:size(edges, 1)
        src = edges(e,1);
        tgt = edges(e,2);
        U{tgt} = U{tgt} + Par(4) * (bs * cs) * errDiff(src, tgt);
    end

    % Error system with reference nodes X_1 and X_2 = X_1 - e12
    Xref = {xs; xs - E{1}};
    F_E = sym([]);
    for p = 1:num_pairs
        i = pairs(p,1);
        j = pairs(p,2);
        F_E = [F_E; Fx(Xref{i}) - Fx(Xref{i} - E{p}) + U{i} - U{j}]; %#ok<AGROW>
    end
    % Independent system: keep de12/dt and de13/dt, and remove e23 from
    % them with the loop closure e23 = e13 - e12 (rows 1:6, every topology).
    % de23/dt is redundant; keeping it adds a fictitious-node mode whose
    % exponent is not physical (e.g. the uncoupled one in the auxiliary motif).
    F_E = subs(F_E(1:dim_V), E{3}, E{2} - E{1});

    % Jacobian w.r.t. [e12; e13] on the synchronization manifold
    J_E0  = subs(jacobian(F_E, es(1:dim_V)), es(1:dim_V), zeros(dim_V, 1));
    F_var = [Fx(xs); J_E0 * dY];

    % Same file and signature used by MaxTransLyaExp
    Zf = str2sym(compose('Z%d(t)', (1:n+dim_V).'));
    odeFunction(subs(F_var, [xs; dY], Zf), Zf, Par, ...
                'File', 'odeMotifVariational.m', 'Optimize', true);
end

% =========================================================================
% LOCAL FUNCTION: EDGE REFINEMENT
%   For every zero crossing of LLE_row (sign change between ks(c) and
%   ks(c+1); NaN counts as not synchronized), the points within
%   +/- half (ks units) are recomputed with opts_ref (longer tf). If, after
%   refining, the whole window has the sign of one side, the crossing has
%   moved out of it and the window is widened by half on that side (at
%   most max_expand times). Points are never computed twice.
% =========================================================================
function [row, done] = refineEdges(row, ks_vec, OscPar, dim_Y, n, opts_ref, half, max_expand, use_parallel)
    N    = numel(row);
    done = false(1, N);
    w    = max(1, round(half / (ks_vec(2) - ks_vec(1))));
    sgn  = row < 0;                                 % Synchronized (NaN -> false)
    cr   = find(diff(sgn) ~= 0);
    for c = cr
        left_val  = sgn(c);                         % Coarse sign on the left of the crossing
        lo = max(1, c - w + 1);
        hi = min(N, c + w);
        lo_all = lo;                                % Whole refined range around this crossing
        hi_all = hi;
        for e = 0:max_expand
            idx = lo:hi;
            idx = idx(~done(idx));
            vals = NaN(1, numel(idx));
            if use_parallel
                parfor kk = 1:numel(idx)
                    vals(kk) = MaxTransLyaExp([OscPar, ks_vec(idx(kk))], dim_Y, n, opts_ref);
                end
            else
                for kk = 1:numel(idx)
                    vals(kk) = MaxTransLyaExp([OscPar, ks_vec(idx(kk))], dim_Y, n, opts_ref);
                end
            end
            row(idx)  = vals;
            done(idx) = true;
            s_all = row(lo_all:hi_all) < 0;
            if any(diff(s_all) ~= 0)
                break;                              % Crossing bracketed by refined points
            elseif all(s_all == ~left_val) && lo_all > 1
                hi = lo_all - 1;                    % Crossing moved to the left: widen on the left
                lo = max(1, lo_all - w);
                lo_all = lo;
            elseif all(s_all == left_val) && hi_all < N
                lo = hi_all + 1;                    % Crossing moved to the right: widen on the right
                hi = min(N, hi_all + w);
                hi_all = hi;
            else
                break;                              % Reached the end of the ks range
            end
        end
    end
end

% =========================================================================
% LOCAL FUNCTION: WORKERS USE THE CLIENT FOLDER AND THE NEW ODE FILE
%   The pool may have been started from another folder (its workers would
%   then find another odeMotifVariational.m on their path), and workers
%   cache the previously loaded version: move every worker to the client's
%   folder, clear the cached function and rescan the path.
% =========================================================================
function refreshWorkers()
    if isempty(gcp('nocreate'))
        return;
    end
    here = strrep(pwd, '''', '''''');
    pctRunOnAll(['cd(''' here ''')']);
    pctRunOnAll('clear odeMotifVariational');
    pctRunOnAll('rehash');
end

% =========================================================================
% LOCAL FUNCTION: ODE FILE, FOLDER AND RANDOM GENERATOR SEEN BY A WORKER
% =========================================================================
function info = workerInfo()
    r    = rng;
    info = sprintf('%s (folder %s, rng %s)', which('odeMotifVariational'), pwd, r.Type);
end

% =========================================================================
% LOCAL FUNCTION: TEXT CHOICE
% =========================================================================
function str = ternaryText(cond, str_true, str_false)
    if cond
        str = str_true;
    else
        str = str_false;
    end
end

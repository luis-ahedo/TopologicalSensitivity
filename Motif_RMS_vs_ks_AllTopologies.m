% =========================================================================
% SCRIPT: Motif_RMS_vs_ks_AllTopologies.m
% DESCRIPTION: Computes rms_e(ks) = mean([rms(e12), rms(e13), rms(e23)])
%              over the steady state (transient removed) for ks in
%              I_ks = [0, kf], for CasoBs = 1 to 9 and the topologies in
%              topology_run (all 13 by default). No figures (they are
%              made by RMS_PlotAndTable.m). Detects the ks intervals where
%              rms_e < sync_tol and prints them. With early_stop, a run
%              stops once the error stays below e_stop for t_chunk time
%              units (the rest of the RMS window counts as zero).
%              Each topology is saved in its own file as soon as it
%              finishes,
%                  mat/RMS_vs_ks_<Osc>_<topology><time_tag>.mat,
%              so the topologies can be computed in groups (e.g. per
%              night). RMS_Consolidate.m joins them into
%              mat/RMS_SyncIntervals_<Osc>_AllTopologies<time_tag>.mat, from
%              which RMS_PlotAndTable.m makes the figures and LaTeX
%              tables; RMS_Consolidate.m runs at the end of this
%              script when consolidate_at_end = true.
% =========================================================================
clear; clc; close all;

% =========================================================================
% STEP 1: INITIALIZATION & OSCILLATOR SETUP
% =========================================================================
OsID = 1;               % 1=Rossler, 2=Lorenz
num_oscillators = 3;

syms a b c ks t real
[F, OscPar, FMtype] = OscillatorType(OsID);

syms x(t) y(t) z(t)
X = [x(t); y(t); z(t)];
Par = [a, b, c, ks];
n = length(X);

% =========================================================================
% STEP 2: SWEEP & SIMULATION CONFIGURATION
% =========================================================================
% ks grid: I_ks = [0, kf] with step dks, kf chosen from the oscillator
%                Rossler  Lorenz  Chen   HR    Duffing  vanderPol
kf_table   = [   5,       25,     25,    25,   25,      25    ];
dks        = 0.01;                          % ks step (two-decimal resolution)
kf         = kf_table(OsID);
N_ks       = round(kf / dks) + 1;           % 501 (Rossler), 2501 (Lorenz)
ks_vec     = linspace(0, kf, N_ks);

% Time horizon profile (same names and values as SensitivityNorm_vs_ks.m,
% so both scripts measure over the same window):
%   'original' : t_transient =  100, t_end =  500   (files without suffix)
%   'short'    : t_transient =  200, t_end = 1000   (files ..._T1000)
%   'medium'   : t_transient =  500, t_end = 2500   (files ..._T2500)
%   'long'     : t_transient = 1000, t_end = 5000   (files ..._T5000)
time_profile = 'medium';
switch lower(time_profile)
    case 'original', t_transient = 100;  t_end = 500;
    case 'short',    t_transient = 200;  t_end = 1000;
    case 'medium',   t_transient = 500;  t_end = 2500;
    case 'long',     t_transient = 1000; t_end = 5000;
    otherwise, error('time_profile must be ''original'', ''short'', ''medium'' or ''long''.');
end
if strcmpi(time_profile, 'original')
    time_tag = '';
else
    time_tag = sprintf('_T%d', t_end);
end
tspan       = 0:0.01:t_end;         % Uniform output grid for the rms (discard t < t_transient)
sync_tol    = 1e-3;                 % rms_e < sync_tol  =>  synchronized

% Convergence stop: the run is integrated in chunks of t_chunk time units;
% once a whole chunk has every ||e_ij|| below e_stop (before or inside the
% RMS window), the run stops and the rest of the window counts as zero
% error. That changes rms_e by less than e_stop, 1000 times below sync_tol.
early_stop = true;
e_stop     = 1e-3 * sync_tol;       % 1e-6: error bound considered "converged"
t_chunk    = 50;                    % Chunk length = time the error must stay below e_stop

% Topologies computed in this run (indices 1..13 of topology_list below),
% e.g. 1:4 tonight, 5:9 tomorrow, 10:13 the next night.
%   1 chain, 2 auxiliar, 3 relay, 4 mixed_master, 5 competitive,
%   6 mixed_slave, 7 ring_bidir_chain, 8 ring_bidir_double_mixed,
%   9 ring_bidir_mixed, 10 ring_unidir_chain, 11 ring_bidir_auxiliar,
%   12 ring_unidir_auxiliar, 13 ring_bidir_competitive
topology_run = 1:13;

% At the end of the run, join the per-topology files with RMS_Consolidate.m
% (it only writes the consolidated file when all 13 topologies exist,
% unless its allow_incomplete = true)
consolidate_at_end = true;

max_state_bound = 1e4;
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, ...
                 'Events', @(t, Y) blowupEvent(t, Y, max_state_bound));


% =========================================================================
% STEP 3: TOPOLOGY DEFINITIONS
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

osc_names = {'R\"{o}ssler', 'Lorenz'};
motif_names = {'Linear chain', 'Linear auxiliary', 'Linear relay', ...
    'Linear mixed master', 'Linear competitive', 'Linear mixed slave', ...
    'Ring bidirectional chain', 'Ring bidirectional double mixed', ...
    'Ring bidirectional mixed', 'Ring unidirectional chain', ...
    'Ring bidirectional auxiliary', 'Ring unidirectional auxiliary', ...
    'Ring bidirectional competitive'};

% =========================================================================
% STEP 4: STATE & VECTOR FIELD GENERATION
% =========================================================================
Xj  = [];
FXj = [];
for k = 1:num_oscillators
    eval(sprintf('syms x%d(t) y%d(t) z%d(t)', k, k, k));
    eval(sprintf('X%d = [x%d(t); y%d(t); z%d(t)];', k, k, k, k));
    eval(sprintf('FX%d = subs(F, X, X%d);', k, k));
    Xj  = [Xj; eval(sprintf('X%d', k))];
    FXj = [FXj; eval(sprintf('FX%d', k))];
end

total_states = n * num_oscillators;
idx = cell(num_oscillators, 1);
for k = 1:num_oscillators
    idx{k} = (k-1)*n + (1:n);
end

rng(42);
initial_conditions = rand(total_states, 1);   % Same ICs for every ks and CasoBs

pairs_matrix = nchoosek(1:num_oscillators, 2);
num_pairs = size(pairs_matrix, 1);
ss_mask = tspan >= t_transient;
N_win   = nnz(ss_mask);                       % samples of the RMS window
% Chunk boundaries as sample indices of tspan (no floating-point gaps)
dt_out = tspan(2) - tspan(1);
if early_stop
    chunk_edges = unique([0:t_chunk:t_end, t_end]);
else
    chunk_edges = [0, t_end];                 % one piece: full integration
end
chunk_idx = round(chunk_edges / dt_out) + 1;  % indices into tspan
i_win     = find(ss_mask, 1);                 % first sample of the RMS window

mat_folder = 'mat';
if ~exist(mat_folder, 'dir')
    mkdir(mat_folder);
end
tex_folder = 'tex';
if ~exist(tex_folder, 'dir')
    mkdir(tex_folder);
end

% =========================================================================
% SCHEME NAMES AND RUN CHECK
% =========================================================================
num_topo = numel(topology_list);
if any(topology_run < 1 | topology_run > num_topo)
    error('topology_run must contain indices between 1 and %d.', num_topo);
end
scheme_names = strings(1, 9);
for CasoBs = 1:9
    [bs, cs] = SISObc(CasoBs, n);
    scheme_names(CasoBs) = sprintf('b%d_to_c%d', find(bs), find(cs));
end

% =========================================================================
% LOOP OVER THE TOPOLOGIES OF THIS RUN
% =========================================================================
t_total = tic;
for topology_idx = topology_run
    t_topo = tic;
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

    fprintf('\n=== Topology %d (%d topologies in this run): %s ===\n', topology_idx, numel(topology_run), motif_names{topology_idx});

    rms_e_all  = NaN(9, N_ks);
    t_stop_all = NaN(9, N_ks);                  % time each run stopped (t_end = full integration)

    % =========================================================================
    % STEP 5: LOOP OVER COUPLING SCHEMES (CasoBs = 1 to 9)
    % =========================================================================
    for CasoBs = 1:9
        [bs, cs] = SISObc(CasoBs, n);
        b_idx = find(bs);
        c_idx = find(cs);
        fprintf('  CasoBs = %d / 9 (b_%d -> c_%d)\n', CasoBs, b_idx, c_idx);

        % Global coupling matrix
        Mc = sym(zeros(total_states, total_states));
        for e = 1:num_edges
            src = edges(e,1);
            tgt = edges(e,2);
            Mc(idx{tgt}, idx{src}) = Mc(idx{tgt}, idx{src}) + ks * (bs * cs);
            Mc(idx{tgt}, idx{tgt}) = Mc(idx{tgt}, idx{tgt}) - ks * (bs * cs);
        end
        F_coupled = FXj + Mc * Xj;

        odeFunction(F_coupled, Xj, Par, 'File', 'odeMotifCoupled.m', 'Optimize', true);
        clear odeMotifCoupled;
        rehash;

        % ---------------------------------------------------------------------
        % Sweep ks in I_ks = [0, kf]
        % ---------------------------------------------------------------------
        for k_idx = 1:N_ks
            Parval = [OscPar(1), OscPar(2), OscPar(3), ks_vec(k_idx)];
            f_ode  = @(t, Y) odeMotifCoupled(t, Y, Parval);

            % Integration in chunks with the convergence stop; the RMS of
            % ||e_ij|| is accumulated over the window samples on the fly
            sumsq   = zeros(1, num_pairs);
            Y0      = initial_conditions;
            blown   = false;
            t_stop  = t_end;
            for s_idx = 1:numel(chunk_edges) - 1
                i_seg = chunk_idx(s_idx):chunk_idx(s_idx + 1);
                t_seg = tspan(i_seg);
                [t_out, Y_out, te] = ode45(f_ode, t_seg, Y0, options);
                if ~isempty(te) || numel(t_out) < numel(t_seg)
                    blown = true;               % Blow-up: leave NaN
                    break;
                end
                Y0 = Y_out(end, :).';
                if s_idx > 1                    % First sample = last of the previous chunk
                    i_seg = i_seg(2:end);
                    Y_out = Y_out(2:end, :);
                end
                in_w  = i_seg(:) >= i_win;
                e_max = 0;
                for p_idx = 1:num_pairs
                    Ei   = Y_out(:, idx{pairs_matrix(p_idx,1)});
                    Ej   = Y_out(:, idx{pairs_matrix(p_idx,2)});
                    e_ij = sqrt(sum((Ei - Ej).^2, 2));  % ||e_ij(t)||
                    sumsq(p_idx) = sumsq(p_idx) + sum(e_ij(in_w).^2);
                    e_max = max(e_max, max(e_ij));
                end
                if early_stop && e_max < e_stop
                    t_stop = chunk_edges(s_idx + 1);  % Converged: rest of the window ~ 0
                    break;
                end
            end
            if blown
                continue;
            end
            t_stop_all(CasoBs, k_idx) = t_stop;
            rms_pairs = sqrt(sumsq / N_win);            % rms over ALL window samples
            rms_e_all(CasoBs, k_idx) = mean(rms_pairs); % (rms12 + rms13 + rms23)/3
        end
        fprintf('    early stop at %d of %d ks values (mean stop time %.0f of %g)\n', ...
            nnz(t_stop_all(CasoBs, :) < t_end), N_ks, mean(t_stop_all(CasoBs, :), 'omitnan'), t_end);
    end


    % ---------------------------------------------------------------------
    % Synchronization intervals: rms_e(ks) < sync_tol
    % ---------------------------------------------------------------------
    sync_intervals = cell(9, 1);        % Each cell: [ks_start, ks_end] rows
    CasoBs_col  = (1:9)';
    scheme_col  = strings(9, 1);
    interval_col = strings(9, 1);
    for CasoBs = 1:9
        [bs, cs] = SISObc(CasoBs, n);
        scheme_col(CasoBs) = sprintf('b%d -> c%d', find(bs), find(cs));
        sync_intervals{CasoBs} = syncIntervals(ks_vec, rms_e_all(CasoBs, :), sync_tol);
        interval_col(CasoBs) = intervalString(sync_intervals{CasoBs}, '%.2f');
    end
    SyncTable = table(CasoBs_col, scheme_col, interval_col, ...
        'VariableNames', {'CasoBs', 'Scheme', 'Sync_intervals_ks'});

    fprintf('\nSynchronization intervals (rms_e < %.0e) - %s, %s:\n', ...
        sync_tol, FMtype, motif_names{topology_idx});
    disp(SyncTable);

    % One file per topology, written as soon as it finishes
    RMS_topo = rms_e_all;                           % 9 x N_ks: RMS_topo(CasoBs, ks index)
    RMS_data = struct();
    for CasoBs = 1:9
        RMS_data.(char(scheme_names(CasoBs))) = rms_e_all(CasoBs, :);
    end
    topology_name = edges_input;
    motif_name    = motif_names{topology_idx};
    run_time_min  = toc(t_topo) / 60;
    date_computed = datestr(now, 'yyyy-mm-dd HH:MM');
    filename_mat  = fullfile(mat_folder, sprintf('RMS_vs_ks_%s_%s%s.mat', FMtype, topology_name, time_tag));
    save(filename_mat, 'RMS_topo', 'RMS_data', 't_stop_all', 'early_stop', 'e_stop', 't_chunk', ...
        'sync_intervals', 'SyncTable', 'interval_col', ...
        'ks_vec', 'sync_tol', 'kf', 'N_ks', 'tspan', 't_transient', 'time_profile', 'time_tag', ...
        'scheme_names', 'topology_idx', 'topology_name', 'motif_name', 'edges', ...
        'osc_names', 'FMtype', 'OsID', 'OscPar', 'run_time_min', 'date_computed');
    fprintf('Topology finished in %.1f min. Saved: %s\n', run_time_min, filename_mat);
end
fprintf('\nTopologies %s finished in %.1f min\n', mat2str(topology_run), toc(t_total)/60);

% =========================================================================
% CONSOLIDATED FILE (RMS_Consolidate.m, same oscillator
% and time profile)
% =========================================================================
if consolidate_at_end
    if exist('RMS_Consolidate.m', 'file')
        fprintf('\nConsolidating the per-topology files...\n');
        consolidate_OsID         = OsID;
        consolidate_time_profile = time_profile;
        run('RMS_Consolidate.m');
    else
        warning('RMS_Consolidate.m not found in the path: consolidated file not written.');
    end
end

% =========================================================================
% LOCAL EVENT FUNCTION FOR ODE45 BLOW-UP TERMINATION
% =========================================================================
function [value, isterminal, direction] = blowupEvent(~, Y, threshold)
    if any(isnan(Y)) || any(isinf(Y))
        value = 0;
    else
        value = max(abs(Y)) - threshold;
    end
    isterminal = 1;
    direction  = 1;
end

% =========================================================================
% LOCAL FUNCTIONS: SYNCHRONIZATION INTERVAL DETECTION
% =========================================================================
function I = syncIntervals(ks_vec, rms_e, tol)
    % Contiguous runs of ks where rms_e < tol (NaN = blow-up, not synced).
    % Returns an N x 2 array [ks_start, ks_end] on the ks grid.
    mask  = rms_e < tol;
    d     = diff([0, mask, 0]);
    i_on  = find(d == 1);
    i_off = find(d == -1) - 1;
    I = [ks_vec(i_on).', ks_vec(i_off).'];
end

function str = intervalString(I, fmt, sep)
    % Formats intervals as "[a, b]; [c, d]" (or with a custom separator).
    if nargin < 3, sep = '; '; end
    if isempty(I)
        str = "--";
        return;
    end
    parts = strings(size(I, 1), 1);
    for r = 1:size(I, 1)
        parts(r) = sprintf(['[' fmt ', ' fmt ']'], I(r, 1), I(r, 2));
    end
    str = join(parts, sep);
end

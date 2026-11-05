% =========================================================================
% SCRIPT: SensitivityIndexMonteCarlo.m
% DESCRIPTION: Self-contained (no other script needs to be run first).
%              Part 1 builds the error system with mismatches and its
%              sensitivity matrix equation symbolically (same code as
%              MotifSensitivityMismatch.m) and generates
%              odeMotifSensitivityMismatch.m for the OsID, topology_idx and
%              CasoBs set in Part 1.
%              Part 2 computes the sensitivity coefficient index (dB),
%              averaged over initial conditions (Monte Carlo), with
%              ALL mismatches = 0 during the integration:
%                 sbar_ij(t) = s_ij(t) * p_j / x_ref,  p_j = nominal a, b or c,
%                                     x_ref = attractor scale (AttractorScale.m)
%                 <sbar_ij>  = sqrt( (1/N) sum_k sbar_ij(t_k)^2 )
%                 S_dB(i,j)  = 20 log10 <sbar_ij>
%              The whole 84-state system is integrated with a fixed-step
%              RK4, several runs at once (one column per run); the step
%              dt is chosen from the oscillator (solver = 'RK4'), or with
%              ode45, one run at a time (solver = 'ode45', RelTol 1e-6,
%              AbsTol 1e-8). Runs stop early once the
%              index converges (at most max_runs), using the mean (SEM) or
%              the median as statistic (stat_predictor); a blow-up test
%              stops any run whose synchronization errors exceed err_bound.
%              Node 1 can start from an attractor pool (use_pool); the
%              errors start on the synchronization manifold, e(0) = 0
%              (ic_mode = 'manifold'), or from node differences of random
%              points (ic_mode = 'random'). Only the
%              54 error rows are analysed (the X1 rows are identically zero).
%              Z = [x (12); vec(S) (72)], x = [X1; e12; e13; e23],
%              S = dx/dq (12 x 6), q = [da2 db2 dc2 da3 db3 dc3].
% =========================================================================
clear; clc; close all;

% #########################################################################
% PART 1: SYMBOLIC GENERATION OF odeMotifSensitivityMismatch.m
%         (same code as MotifSensitivityMismatch.m; OsID, topology_idx and
%          CasoBs are set in STEP 1 and STEP 2 below)
% #########################################################################
% =========================================================================
% STEP 1: INITIALIZATION & OSCILLATOR SETUP
% =========================================================================
OsID = 1;               % 1=Rossler, 2=Lorenz
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
topology_idx = 1;
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
CasoBs = 1;
% ks is a runtime parameter of the generated ODE file (no ks_val here)
[bs, cs] = SISObc(CasoBs, n);         % Inner coupling vectors
% =========================================================================
% STEP 3: AUTOMATED STATE & VECTOR FIELD GENERATION
% =========================================================================
fprintf('Generating states and individual vector fields for N = %d oscillators...\n', num_oscillators);
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
fprintf('Building the error system de12/dt, de13/dt, de23/dt...\n');

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

% Consistency check (N = 3): de23/dt = de13/dt - de12/dt
% once X_2 = X_1 - e12 and e23 = e13 - e12 are imposed
if num_oscillators == 3
    check = F_matrix_E(7:9) - (F_matrix_E(4:6) - F_matrix_E(1:3));
    check = subs(check, [Xj(idx{2}); E_array{3}], ...
                 [Xj(idx{1}) - E_array{1}; E_array{2} - E_array{1}]);
    if all(isAlways(simplify(check) == 0))
        fprintf('Loop-closure check passed: de23/dt = de13/dt - de12/dt\n');
    else
        warning('Loop-closure check failed.');
    end
end

fprintf('\nError system (%s, %s, CasoBs = %d):\n', FMtype, edges_input, CasoBs);
for p = 1:num_pairs
    fprintf('d%s/dt =\n', pair_names{p});
    disp(F_matrix_E((p-1)*n + (1:n)));
end

% =========================================================================
% SUBSTITUTION OF e23 IN de13/dt:  e23 = e13 - e12
% =========================================================================
F_matrix_E(4:6) = expand(subs(F_matrix_E(4:6), E_array{3}, E_array{2} - E_array{1}));
for r = 4:6
    F_matrix_E(r) = collect(F_matrix_E(r), Par);
end
fprintf('\nde13/dt after substituting e23 = e13 - e12:\n');
disp(F_matrix_E(4:6));

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
fprintf('\nBuilding the sensitivity matrix equation dS/dt = F_x S + F_q...\n');
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

% =========================================================================
% STATE AND PARAMETER ORDERING (for the later simulation)
% =========================================================================
state_names = {'x1', 'y1', 'z1'};
for p = 1:num_pairs
    state_names = [state_names, {[pair_names{p} 'x'], [pair_names{p} 'y'], [pair_names{p} 'z']}]; %#ok<AGROW>
end
param_names = {'da2','db2','dc2','da3','db3','dc3'};
fprintf('\nState vector x (Z(1:%d)): %s\n', Nx, strjoin(state_names, ', '));
fprintf('Sensitivity parameters q: %s\n', strjoin(param_names, ', '));
fprintf('S = reshape(Z(%d:%d), %d, %d); column j = dx/dq_j (q_1 = %s, ..., q_%d = %s)\n', ...
    Nx+1, Nz, Nx, Nq, param_names{1}, Nq, param_names{end});
fprintf('Initial conditions: e23(0) = e13(0) - e12(0), S(0) = 0.\n');
fprintf('Parval order: [da2 db2 dc2 da3 db3 dc3 ks]  (a, b, c fixed at OscPar)\n');

% #########################################################################
% PART 2: MONTE CARLO SENSITIVITY INDEX (vectorized fixed-step RK4)
% #########################################################################

% =========================================================================
% MC STEP 1: CONFIGURATION
% =========================================================================
ks_val      = 2.00;             % Coupling strength
T_end       = 1000;             % Integration time of each run
t_transient = 200;              % RMS window: t >= t_transient
dt_sample   = 0.05;             % Sampling step of the RMS (N samples in the window)

% RK4 step chosen automatically from the oscillator (OsID of Part 1):
%              Rossler  Lorenz  Chen    HR     Duffing  vanderPol
dt_table   = [ 0.01,    0.005,  0.002,  0.01,  0.01,    0.005 ];
dt         = dt_table(OsID);

% --- Integrator of the 84-state system ---
%   'RK4'  : fixed-step RK4 with step dt, all runs of a batch at once
%   'ode45': adaptive ode45, one run at a time, output every dt_sample,
%            tolerances ode_RelTol / ode_AbsTol, blow-up test as an event
solver     = 'ode45';
ode_RelTol = 1e-6;
ode_AbsTol = 1e-8;

max_runs   = 100;               % Maximum Monte Carlo runs (random initial conditions)
batch_size = 10;                % Runs integrated together (one column each)
base_seed  = 42;                % Run r uses rng(base_seed + r)

% --- Statistic used as index and in the convergence test ---
%   'SEM'    : index = MeanRMS_dB = 20 log10(mean of the per-run RMS);
%              uncertainty = 20 log10(1 + SEM/MeanRMS)
%   'Median' : index = Median_dB = median of the per-run 20 log10(RMS);
%              uncertainty = 1.2533 * std(dB values) / sqrt(n)
%              (standard error of the median; robust to a few large runs)
stat_predictor = 'Median';

% --- Initial condition of node 1 ---
%   use_pool = true : node 1 starts on its attractor, taken from a pool of
%                     points of one long run of the uncoupled nominal
%                     oscillator (run r uses pool point r)
%   use_pool = false: node 1 starts at a random point in [0,1]^3
%   The random numbers of each run are drawn with the same seed in both
%   cases, so the two options differ ONLY in X1(0) and can be compared
%   run by run.
use_pool         = true;
pool_T_transient = 500;         % Discarded transient of the pool run
pool_spacing     = 50;          % Time between stored pool points (decorrelation)

% --- Initial condition of the synchronization errors ---
%   ic_mode = 'manifold': e12(0) = e13(0) = e23(0) = 0 (start on the
%             synchronization manifold). With all mismatches = 0 the errors
%             stay exactly 0, so the run-to-run variability comes only from
%             X1(0). The Jacobian blocks that couple S_X1 into s(e) are
%             proportional to e, so S_X1 does not enter s(e); its rows are
%             held at 0 (keep mask), which avoids an overflow of S_X1 in
%             topologies where node 1 is driven.
%   ic_mode = 'random'  : errors are node differences of random points in
%             [0,1]^3 (e23 = e13 - e12), as in the previous version.
ic_mode = 'manifold';

% --- Blow-up test on the synchronization errors e12, e13, e23 ---
% A run is stopped (and excluded) when any |e_ij| component exceeds
% err_bound, or when any entry of its state becomes NaN/Inf.
err_bound   = 1e3;
check_every = 100;              % RK4 steps between blow-up checks (ode45: event, every step)

% --- Convergence test (checked after each batch) ---
min_runs     = 10;              % Valid runs before the test is applied
conv_window  = 10;              % Compare with the index conv_window runs earlier
conv_tol_dB  = 1;               % max |change of the index| over the window (dB)
conv_sem_dB  = 1;               % max uncertainty of the index (dB)
stop_on_convergence = true;     % false: all max_runs are done; the first run count
                                % meeting the test is only reported and marked in the plot

% --- Worst-run check (after the Monte Carlo loop) ---
%   true: the run with the largest excess over the median (in any
%   coefficient) and a typical run (closest to the median) are integrated
%   again with the same seed and solver, and their time series are plotted:
%   a physical burst is smooth and lasts several time units, a numerical
%   spike lasts only a few steps.
plot_worst_run = false;

% --- Coefficients treated as zero (structural or first-order zeros) ---
zero_rel_tol = 1e-6;            % MeanRMS <= zero_rel_tol * max(MeanRMS)  (about -120 dB)

F_nom = F;                      % Symbolic node vector field (from Part 1)
a = OscPar(1);
b = OscPar(2);
c = OscPar(3);
% Sensitivities are evaluated at the nominal point: ALL mismatches = 0
Parval = [zeros(1, 6), ks_val];   % [da2 db2 dc2 da3 db3 dc3 ks], order of the ODE file

font_title    = 18;
font_xlabel   = 20;
font_ylabel   = 20;
font_ticks    = 15;

% --- Figure labels for the paper ---
%   show_title  = false: heat map and convergence figure without title
%                 (the caption of the paper describes the case)
show_title  = true;
%   heat_alpha: opacity of the heat-map colours over a white background
%               (1 = full colours, smaller = softer). The colours are blended
%               with white in the colormap itself, so the EPS stays vector
%               (true transparency would be rasterized or lost in the EPS)
heat_alpha  = 0.55;

osc_names = {'R\"{o}ssler', 'Lorenz'};

% =========================================================================
% MC STEP 2: LABELS AND SCALE VECTOR (order of vec(S): column j = q_j,
%         12 states per column)
% =========================================================================
Nx = 12;
Nq = 6;
Nz = Nx + Nx*Nq;                                  % 84 states
paramNames = {'da2','db2','dc2','da3','db3','dc3'};
stateNames = {'x1','y1','z1','ex12','ey12','ez12','ex13','ey13','ez13','ex23','ey23','ez23'};
err_idx    = 4:12;                                % error rows of x
Ncoef      = numel(err_idx) * Nq;                 % 54 coefficients

% Scale of sbar_ij = s_ij * p_j / x_ref: nominal parameter of each mismatch
% (a for da2 and da3, b for db2 and db3, c for dc2 and dc3)
p_scale = [a b c a b c];
% Attractor scale x_ref of Eq. (eq_normalized_sensitivity), one value per
% oscillator (RMS of the nominal node trajectory, AttractorScale.m):
% sbar_ij = s_ij * q_j / x_ref
xref = AttractorScale(OsID);

% Rows of Z holding the 54 analysed sensitivities (q_j-major) and their scale
sens_rows = zeros(Ncoef, 1);
sens_scal = zeros(Ncoef, 1);
S_labels  = cell(Ncoef, 1);
kk = 1;
for j = 1:Nq
    for i = err_idx
        sens_rows(kk) = Nx + (j-1)*Nx + i;
        sens_scal(kk) = p_scale(j) / xref;          % sbar_ij = s_ij * q_j / x_ref
        S_labels{kk}  = sprintf('S_%s_%s', paramNames{j}, stateNames{i});
        kk = kk + 1;
    end
end

% =========================================================================
% MC STEP 2b: POOL OF NODE-1 INITIAL CONDITIONS ON THE ATTRACTOR (use_pool)
%   One long run of the uncoupled nominal oscillator; after
%   pool_T_transient one point is stored every pool_spacing time units.
% =========================================================================
if use_pool
    xs_pool = sym('xs_pool', [n 1], 'real');
    X_old   = str2sym({'x(t)'; 'y(t)'; 'z(t)'});
    P_old   = [sym('a', 'real'), sym('b', 'real'), sym('c', 'real')];
    F_pool  = subs(F_nom, [X_old; P_old.'], [xs_pool; OscPar(1:3).']);
    f_pool  = matlabFunction(F_pool, 'Vars', {xs_pool});

    t_pool  = [0, pool_T_transient + (1:max_runs) * pool_spacing];
    rng(base_seed);
    [~, Xp] = ode45(@(t, x) f_pool(x), t_pool, rand(n, 1), ...
                    odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
    X1_pool = Xp(2:end, :);                        % max_runs x n, one row per run
    fprintf('Pool of %d node-1 points on the attractor (spacing %g, after t = %g).\n', ...
        max_runs, pool_spacing, pool_T_transient);
else
    X1_pool = NaN(max_runs, n);                    % Not used
end

% =========================================================================
% MC STEP 2c: VECTORIZED RIGHT-HAND SIDE
%   The generated file evaluates several initial conditions at once when
%   they are passed as columns. This is tested once; if the generated code
%   does not support it, the columns are evaluated one by one (same
%   results, slower).
% =========================================================================
% Keep mask: rows of dS_X1/dt held at 0 in ic_mode = 'manifold'
keep = ones(Nz, 1);
if strcmpi(ic_mode, 'manifold')
    for j = 1:Nq
        keep(Nx + (j-1)*Nx + (1:n)) = 0;           % S_X1 rows of column q_j
    end
elseif ~strcmpi(ic_mode, 'random')
    error('ic_mode must be ''manifold'' or ''random''.');
end
f_one = @(Z) keep .* odeMotifSensitivityMismatch(0, Z, Parval);
Ztest = randn(Nz, 3);
vec_ok = false;
try
    Ftest  = f_one(Ztest);
    Fcol   = [f_one(Ztest(:,1)), f_one(Ztest(:,2)), f_one(Ztest(:,3))];
    vec_ok = isequal(size(Ftest), [Nz 3]) && max(abs(Ftest(:) - Fcol(:))) < 1e-10 * max(1, max(abs(Fcol(:))));
catch
    vec_ok = false;
end
if vec_ok
    f_rhs = f_one;
    fprintf('Vectorized evaluation: ON (all runs of a batch in one call).\n');
else
    f_rhs = @(Z) evalColumns(f_one, Z);
    fprintf('Vectorized evaluation: OFF (column by column).\n');
end

% =========================================================================
% MC STEP 3: MONTE CARLO RUNS (BATCHES) WITH CONVERGENCE TEST
%   Each run: node 1 from the pool (or random), errors e(0) = 0
%   (ic_mode = 'manifold') or node differences of random points
%   (e23(0) = e13(0) - e12(0), ic_mode = 'random'), S(0) = 0.
%   RK4 with step dt; the RMS of s_ij * p_j / x_ref is accumulated on the fly
%   every dt_sample for t >= t_transient. After each batch the chosen
%   index (stat_predictor) is compared with its value conv_window runs
%   earlier; the loop stops when every active coefficient changed less than
%   conv_tol_dB and its uncertainty is below conv_sem_dB.
% =========================================================================
N_steps   = round(T_end / dt);
k_trans   = round(t_transient / dt);
k_sample  = max(1, round(dt_sample / dt));
k_q1_end  = k_trans + round((N_steps - k_trans) / 4);     % first quarter of the window
k_q4_beg  = N_steps - round((N_steps - k_trans) / 4);     % last quarter of the window

RMS_runs   = NaN(max_runs, Ncoef);             % per-run RMS (q_j-major)
growth_run = NaN(max_runs, 1);                 % late / early RMS ratio
blown_run  = false(max_runs, 1);               % run stopped by the blow-up test

hist_r     = [];                               % runs completed at each check
hist_idx   = [];                               % index (dB) at each check (rows)
delta_hist = [];                               % max |change| over the window (dB)
semdB_hist = [];                               % max uncertainty (dB)

if strcmpi(solver, 'ode45')
    solver_info = sprintf('ode45, RelTol = %g, AbsTol = %g', ode_RelTol, ode_AbsTol);
    t_samp   = t_transient + (1:floor((T_end - t_transient) / dt_sample)) * dt_sample;
    t_q1_end = k_q1_end * dt;                      % same quarters as RK4 (time units)
    t_q4_beg = k_q4_beg * dt;
    in_q1    = t_samp <= t_q1_end;
    in_q4    = t_samp >= t_q4_beg;
    ode_opts = odeset('RelTol', ode_RelTol, 'AbsTol', ode_AbsTol, ...
                      'Events', @(t, Z) blowupEvent(Z, err_idx, err_bound));
    f_ode    = @(t, Z) f_one(Z);
elseif strcmpi(solver, 'RK4')
    solver_info = sprintf('RK4, dt = %g', dt);
    f_ode    = [];                                 % Not used with RK4
    ode_opts = [];
else
    error('solver must be ''RK4'' or ''ode45''.');
end
fprintf('Monte Carlo sensitivity index: up to %d runs (%s, %s, CasoBs = %d, ks = %.3f, %s, %s, ic_mode = %s)\n', ...
    max_runs, FMtype, topology_list{topology_idx}, CasoBs, ks_val, solver_info, stat_predictor, ic_mode);
t_start   = tic;
r_done    = 0;
converged = false;
r_conv    = [];                                % first run count meeting the test
while r_done < max_runs && ~(converged && stop_on_convergence)
    batch = r_done + 1 : min(r_done + batch_size, max_runs);
    nb    = numel(batch);

    % Initial conditions of the batch (one column per run)
    Zc = zeros(Nz, nb);
    for bb = 1:nb
        rng(base_seed + batch(bb));
        X0 = reshape(rand(n * num_oscillators, 1), n, num_oscillators);
        E0 = [X0(:,1) - X0(:,2); X0(:,1) - X0(:,3); X0(:,2) - X0(:,3)];
        if strcmpi(ic_mode, 'manifold')
            E0 = zeros(size(E0));                  % Start on the synchronization manifold
        end
        if use_pool
            x1_0 = X1_pool(batch(bb), :).';
        else
            x1_0 = X0(:,1);
        end
        Zc(1:Nx, bb) = [x1_0; E0];                 % S(0) = 0
    end

    alive  = true(1, nb);
    sumsq  = zeros(Ncoef, nb);  cnt  = 0;          % whole window
    sumsq1 = zeros(1, nb);      cnt1 = 0;          % first quarter (growth)
    sumsq4 = zeros(1, nb);      cnt4 = 0;          % last quarter  (growth)

    if strcmpi(solver, 'ode45')
    % --- ode45: one run at a time, output only at the RMS sample times ---
    cnt  = numel(t_samp);
    cnt1 = nnz(in_q1);
    cnt4 = nnz(in_q4);
    for bb = 1:nb
        [t_out, Z_out, te] = ode45(f_ode, [0, t_samp], Zc(:, bb), ode_opts);
        if ~isempty(te) || numel(t_out) < cnt + 1 || any(~isfinite(Z_out(:)))
            alive(bb) = false;                     % blow-up event: excluded later
            continue;
        end
        Zs = Z_out(2:end, :).';                    % Nz x samples (t >= t_transient)
        Sb = Zs(sens_rows, :) .* sens_scal;        % sbar_ij, 54 x samples
        sumsq(:, bb) = sum(Sb.^2, 2);
        sumsq1(bb)   = sum(mean(Sb(:, in_q1).^2, 1));
        sumsq4(bb)   = sum(mean(Sb(:, in_q4).^2, 1));
    end
    else
    % --- RK4: all runs of the batch at once ---
    for k = 1:N_steps
        k1 = f_rhs(Zc);
        k2 = f_rhs(Zc + 0.5*dt*k1);
        k3 = f_rhs(Zc + 0.5*dt*k2);
        k4 = f_rhs(Zc + dt*k3);
        Zc = Zc + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);

        % Blow-up test on the errors (and NaN/Inf anywhere)
        if mod(k, check_every) == 0
            bad = any(~isfinite(Zc), 1) | any(abs(Zc(err_idx, :)) > err_bound, 1);
            if any(bad & alive)
                alive(bad) = false;
                Zc(:, bad) = 0;                    % frozen, excluded later
            end
        end

        % On-the-fly RMS accumulation in the window
        if k > k_trans && mod(k - k_trans, k_sample) == 0
            Sb = Zc(sens_rows, :) .* sens_scal;    % sbar_ij, 54 x nb
            sumsq = sumsq + Sb.^2;  cnt = cnt + 1;
            if k <= k_q1_end
                sumsq1 = sumsq1 + mean(Sb.^2, 1);  cnt1 = cnt1 + 1;
            end
            if k >= k_q4_beg
                sumsq4 = sumsq4 + mean(Sb.^2, 1);  cnt4 = cnt4 + 1;
            end
        end
    end
    end

    % Store the batch
    rms_b = sqrt(sumsq / max(cnt, 1));
    rms_b(:, ~alive) = NaN;
    RMS_runs(batch, :)  = rms_b.';
    growth_run(batch)   = sqrt(sumsq4 / max(cnt4, 1)) ./ max(sqrt(sumsq1 / max(cnt1, 1)), eps);
    growth_run(batch(~alive)) = Inf;
    blown_run(batch)    = ~alive.';
    r_done = batch(end);

    % Running statistics over the valid runs so far
    Rv      = RMS_runs(1:r_done, :);
    Rv      = Rv(~blown_run(1:r_done), :);
    n_valid = size(Rv, 1);
    msg = sprintf('  runs %d/%d: %d valid, %d blown up', r_done, max_runs, n_valid, nnz(blown_run(1:r_done)));
    if n_valid > 0
        [idx_now, unc_now] = indexStat(Rv, stat_predictor);
        hist_r   = [hist_r; r_done];                       %#ok<AGROW>
        hist_idx = [hist_idx; idx_now];                    %#ok<AGROW>
        m_now    = mean(Rv, 1);
        active   = m_now > zero_rel_tol * max(m_now);

        % Convergence test against the check conv_window runs earlier
        i_ref = find(hist_r <= r_done - conv_window, 1, 'last');
        if n_valid >= min_runs && ~isempty(i_ref)
            delta  = max(abs(idx_now(active) - hist_idx(i_ref, active)));
            sem_dB = max(unc_now(active));
            delta_hist = [delta_hist; r_done, delta];      %#ok<AGROW>
            semdB_hist = [semdB_hist; r_done, sem_dB];     %#ok<AGROW>
            msg = sprintf('%s | max dIndex = %.3f dB, max uncertainty = %.3f dB', msg, delta, sem_dB);
            if delta < conv_tol_dB && sem_dB < conv_sem_dB
                msg = [msg ' -> test met'];                %#ok<AGROW>
                if isempty(r_conv)
                    r_conv = r_done;
                end
            end
            converged = ~isempty(r_conv);
        end
    end
    fprintf('%s  (%.1f min)\n', msg, toc(t_start)/60);
end
fprintf('Finished in %.1f min after %d runs', toc(t_start)/60, r_done);
if converged
    fprintf(' (converged at %d runs, %s).\n', r_conv, stat_predictor);
else
    fprintf(' (max_runs reached, NOT converged).\n');
end

% Keep only the runs actually performed
RMS_runs   = RMS_runs(1:r_done, :);
growth_run = growth_run(1:r_done);
blown_run  = blown_run(1:r_done);
if nnz(blown_run) > 0
    warning('%d of %d runs blew up (|e_ij| > %g) and were excluded.', nnz(blown_run), r_done, err_bound);
end

% =========================================================================
% MC STEP 4: STATISTICS OVER THE VALID RUNS
%   MeanRMS_dB : 20 log10 of the mean linear RMS
%   Median_dB, P05_dB, P95_dB : percentiles of the per-run dB values
%   Index_dB   : the one selected by stat_predictor
% =========================================================================
R2      = RMS_runs(~blown_run, :);             % valid runs x 54 (q_j-major)
n_valid = size(R2, 1);
S_mean  = mean(R2, 1, 'omitnan');
S_std   = std(R2, 0, 1, 'omitnan');
S_sem   = S_std ./ sqrt(sum(~isnan(R2), 1));   % Standard error of the mean
dB_runs = 20*log10(R2);
[Index_dB, Unc_dB] = indexStat(R2, stat_predictor);
S_table = table(S_labels, S_mean(:), S_std(:), S_sem(:), 20*log10(S_mean(:)), ...
    median(dB_runs, 1, 'omitnan').', colPrctile(dB_runs, 5).', colPrctile(dB_runs, 95).', ...
    Index_dB(:), Unc_dB(:), ...
    'VariableNames', {'Quantity', 'MeanRMS', 'StdRMS', 'SEM', 'MeanRMS_dB', ...
                      'Median_dB', 'P05_dB', 'P95_dB', 'Index_dB', 'Uncert_dB'});

% Remove the zero coefficients (below zero_rel_tol * max, about -120 dB):
% - structural zeros: the mismatch has no path to that error in this
%   topology, e.g. S_da3_ex12 in the linear chain (e12 does not depend on node 3)
% - first-order zeros: e.g. S_da2_e23 for Lorenz with x-coupling, where the
%   sigma mismatch only time-shifts node 2 and node 3 follows it exactly
zero_coef = S_mean <= zero_rel_tol * max(S_mean);
if any(zero_coef)
    fprintf('\nRemoved %d zero coefficients (< %g x max): %s\n', nnz(zero_coef), zero_rel_tol, ...
        strjoin(S_labels(zero_coef).', ', '));
end
S_table = S_table(~zero_coef, :);

fprintf('\nSensitivity indices over %d valid runs (scaled by nominal a, b, c; Index_dB = %s):\n', ...
    n_valid, stat_predictor);
disp(S_table);

% Growth check: if the sensitivities keep growing (synchronization not
% stable), the index depends on the window length and must be read with care
growth_med = median(growth_run(~blown_run), 'omitnan');   % valid runs only
fprintf('Median late/early RMS ratio over valid runs: %.2f', growth_med);
if growth_med > 10
    fprintf('  -> sensitivities are GROWING: the index depends on T_end.\n');
else
    fprintf('  -> sensitivities are bounded over the window.\n');
end

% =========================================================================
% MC STEP 4b: RMS SENSITIVITY VECTOR INDEX OF EACH ERROR VECTOR
%   s_qj(e_ij)   = [s_qj(e_x), s_qj(e_y), s_qj(e_z)]^T   (linear units)
%   ||s_qj(e_ij)|| = ( s_qj(e_x)^2 + s_qj(e_y)^2 + s_qj(e_z)^2 )^(1/2)
%   Components: the index of MC STEP 4 in linear units (median or mean of
%   the per-run RMS of sbar_ij, as stat_predictor); zero coefficients = 0.
%   Norm uncertainty: the SAME rule as Uncert_dB of the components
%   (indexStat), applied to the per-run norms: 'Median' -> 1.2533 x std of
%   the per-run norms in dB / sqrt(n); 'SEM' -> 20 log10(1 + SEM/mean).
%   V_unc_dB holds it in dB; V_std is the same uncertainty in linear units.
% =========================================================================
idx_lin = 10.^(Index_dB / 20);                     % 1 x 54 (q_j-major)
idx_lin(zero_coef) = 0;
pair_lbl = {'e12', 'e13', 'e23'};
nV       = Nq * num_pairs;                         % 18 sensitivity vectors
V_label  = cell(nV, 1);
V_param  = cell(nV, 1);
V_err    = cell(nV, 1);
V_comp   = zeros(nV, n);
V_norm   = zeros(nV, 1);
V_std    = zeros(nV, 1);                          % linear uncertainty of the norm
V_unc_dB = zeros(nV, 1);                          % same uncertainty in dB
spread_name = 'uncertainty';
row = 0;
for j = 1:Nq
    for p = 1:num_pairs
        row  = row + 1;
        cols = (j-1)*numel(err_idx) + (p-1)*n + (1:n);        % [e_x e_y e_z] of pair p
        V_comp(row, :) = idx_lin(cols);
        V_norm(row)    = norm(idx_lin(cols));
        if ~all(zero_coef(cols))
            nrm_runs = sqrt(sum(R2(:, cols).^2, 2));              % per-run norm (linear)
            [~, V_unc_dB(row)] = indexStat(nrm_runs, stat_predictor);  % same rule as the components
            V_std(row) = V_norm(row) * (10^(V_unc_dB(row) / 20) - 1);
        end
        V_param{row} = paramNames{j}(2:end);                  % 'a2', 'b2', ...
        V_err{row}   = pair_lbl{p};
        V_label{row} = sprintf('s_%s(%s)', V_param{row}, V_err{row});
    end
end
V_table = table(V_label, V_param, V_err, V_comp(:,1), V_comp(:,2), V_comp(:,3), V_norm, V_std, ...
    'VariableNames', {'SensitivityVector', 'Parameter', 'ErrorVector', 'e_x', 'e_y', 'e_z', ...
                      'Norm', 'Norm_Uncert'});
% First table in dB (the index is defined in dB): components = Index_dB,
% norm = 20 log10 of the norm computed on the LINEAR components, and its
% uncertainty in dB from indexStat (same rule as Uncert_dB). Exact zeros = -Inf.
V_comp_dB = 20 * log10(V_comp);
V_norm_dB = 20 * log10(V_norm);
V_table_dB = table(V_label, V_param, V_err, V_comp_dB(:,1), V_comp_dB(:,2), V_comp_dB(:,3), V_norm_dB, V_unc_dB, ...
    'VariableNames', {'SensitivityVector', 'Parameter', 'ErrorVector', 'e_x_dB', 'e_y_dB', 'e_z_dB', ...
                      'Norm_dB', 'Norm_Uncert_dB'});
fprintf('\nRMS sensitivity vector index of each error vector (%s over %d runs, dB):\n', ...
    stat_predictor, n_valid);
fprintf('%s\n', repmat('-', 1, 110));
fprintf('%-20s %-10s %-13s %10s %10s %10s %12s %20s\n', 'Sensitivity vector', 'Parameter', ...
    'Error vector', 'e_x', 'e_y', 'e_z', '||s_qj(e)||', ['+/- ' spread_name]);
fprintf('%s\n', repmat('-', 1, 110));
for row = 1:nV
    fprintf('%-20s %-10s %-13s %10s %10s %10s %12s %20s\n', V_label{row}, V_param{row}, V_err{row}, ...
        dbStr(V_comp_dB(row, 1), false), dbStr(V_comp_dB(row, 2), false), dbStr(V_comp_dB(row, 3), false), ...
        dbStr(V_norm_dB(row), false), dbStr(V_unc_dB(row), false));
    if mod(row, num_pairs) == 0 && row < nV
        fprintf('\n');                                   % Blank line between parameters
    end
end
fprintf('%s\n', repmat('-', 1, 110));

% Alternative layout (kept for comparison): same as the heat map.
% Same layout as the heat map: columns = mismatches q_j, rows = error
% vectors e12, e13, e23 with their components x, y, z and the norm
% (+/- its uncertainty) below each vector.
comp_lbl = {'x', 'y', 'z', '||.||'};
V_grid     = zeros(num_pairs * (n + 1), Nq);          % 12 x 6: [x; y; z; norm] per error vector
V_grid_unc = NaN(num_pairs * (n + 1), Nq);            % uncertainty on the norm rows only
for p = 1:num_pairs
    for j = 1:Nq
        rv = (j - 1) * num_pairs + p;                  % row of V_comp / V_norm / V_std
        rg = (p - 1) * (n + 1);
        V_grid(rg + (1:n), j)  = V_comp(rv, :).';
        V_grid(rg + n + 1, j)  = V_norm(rv);
        V_grid_unc(rg + n + 1, j) = V_std(rv);
    end
end
fmt_norm = @(v, u) ternary_str(v == 0, '0', sprintf('%.4g +/- %.2g', v, u));
w_col = 18;
line_w = 14 + w_col * Nq;
fprintf('\nSame data, heat-map layout (%s over %d runs, linear units):\n', ...
    stat_predictor, n_valid);
fprintf('%s\n', repmat('-', 1, line_w));
fprintf('%-7s %-6s', 'Error', 'Comp');
for j = 1:Nq
    fprintf('%*s', w_col, paramNames{j});
end
fprintf('\n%s\n', repmat('-', 1, line_w));
for p = 1:num_pairs
    for c = 1:n + 1
        rg = (p - 1) * (n + 1) + c;
        if c == 1, e_lbl = pair_lbl{p}; else, e_lbl = ''; end
        fprintf('%-7s %-6s', e_lbl, comp_lbl{c});
        for j = 1:Nq
            if c <= n
                fprintf('%*s', w_col, sprintf('%.4g', V_grid(rg, j)));
            else
                fprintf('%*s', w_col, fmt_norm(V_grid(rg, j), V_grid_unc(rg, j)));
            end
        end
        fprintf('\n');
    end
    if p < num_pairs, fprintf('\n'); end
end
fprintf('%s\n', repmat('-', 1, line_w));

% =========================================================================
% MC STEP 5: HEAT MAP OF THE INDEX (9 error components x 6 mismatches)
% =========================================================================
S_dB_map = reshape(Index_dB, numel(err_idx), Nq);   % rows: error comps
comp_tex = {'e_{x_{12}}','e_{y_{12}}','e_{z_{12}}','e_{x_{13}}','e_{y_{13}}','e_{z_{13}}', ...
            'e_{x_{23}}','e_{y_{23}}','e_{z_{23}}'};      % Notation of the paper
par_sym = {'a', 'b', 'c'};                          % Plain letters (MATLAB has no lowercase \mathcal)
par_tex = cell(1, Nq);
for j = 1:Nq
    par_tex{j} = sprintf('\\delta %s_%d', par_sym{mod(j-1, 3) + 1}, 2 + (j > 3));   % \delta a_2 ... \delta c_3
end
[bs, cs] = SISObc(CasoBs, n);

fig = figure('Name', 'Sensitivity index (dB)');
ax  = axes(fig);
zero_map = reshape(zero_coef, numel(err_idx), Nq);  % Structural/first-order zeros
% Zero cells in grey: a grey entry is added below the colormap and the zero
% cells get a value that maps to it, so the grey is part of the image itself
% (kept in the EPS). The colour limits come from the non-zero coefficients.
cmap_n   = 256;
cmap     = nebula(cmap_n);
cmap     = heat_alpha * cmap + (1 - heat_alpha);     % Softer colours: blend with white
act_vals = S_dB_map(~zero_map);
if isempty(act_vals)
    c_lo = 0; c_hi = 1;
else
    c_lo = floor(min(act_vals));
    c_hi = max(ceil(max(act_vals)), c_lo + 1);
end
c_step   = (c_hi - c_lo) / cmap_n;                  % Width of one colour
map_plot = S_dB_map;
map_plot(zero_map) = c_lo - c_step / 2;             % -> grey entry
imagesc(ax, map_plot);
colormap(ax, [0.85 0.85 0.85; cmap]);
set(ax, 'CLim', [c_lo - c_step, c_hi]);
cb = colorbar(ax);
cb.Limits = [c_lo, c_hi];                           % Grey entry hidden in the colorbar
cb.TickLabelInterpreter = 'latex';
cb.Label.Interpreter = 'latex';
if strcmpi(stat_predictor, 'Median')
    cb.Label.String = '$s_{q_j}(x_i) = \mathrm{median}_r\, 20\log_{10}\mathrm{rms}(\bar{s}_{i,j}^{(r)})$ (dB)';
else
    cb.Label.String = '$20\log_{10}[\mathrm{mean}_r\, \mathrm{rms}(\bar{s}_{i,j}^{(r)})]$ (dB)';
end
cb.Label.FontSize = font_ylabel;
set(ax, 'FontSize', font_ticks, 'TickLabelInterpreter', 'latex', ...
    'XTick', 1:Nq, 'XTickLabel', cellfun(@(s) sprintf('$%s$', s), par_tex, 'UniformOutput', false), ...
    'YTick', 1:numel(err_idx), 'YTickLabel', cellfun(@(s) sprintf('$%s$', s), comp_tex, 'UniformOutput', false));
clim_ax = get(ax, 'CLim');                        % Text colour readable on the colormap
for i = 1:numel(err_idx)                     % Values printed in each cell
    for j = 1:Nq
        if zero_map(i, j)
            txt_col = 'k';                   % Zero coefficient (grey cell)
            txt_str = '$-\infty$';                % dB of an exactly zero coefficient
        else
            % Text colour from the luminance of the (blended) cell colour
            k_c = round((S_dB_map(i, j) - c_lo) / (c_hi - c_lo) * (cmap_n - 1)) + 1;
            k_c = min(max(k_c, 1), cmap_n);
            lum = cmap(k_c, :) * [0.299; 0.587; 0.114];
            if lum > 0.5
                txt_col = 'k';               % Light cells
            else
                txt_col = 'w';               % Dark cells
            end
            txt_str = sprintf('$%.2f$', S_dB_map(i, j));
        end
        txt_str = strrep(txt_str, '$-0.00$', '$0.00$');   % No '-0.00' for tiny negative values
        text(ax, j, i, txt_str, 'HorizontalAlignment', 'center', 'Interpreter', 'latex', ...
            'FontSize', font_ticks - 3, 'Color', txt_col);
    end
end
xlabel(ax, 'Mismatch $q_j$', 'Interpreter', 'latex', 'FontSize', font_xlabel);
ylabel(ax, 'Error component $x_i$', 'Interpreter', 'latex', 'FontSize', font_ylabel);
if use_pool, ic_str = ', $\mathbf{x}_1(0)$ from pool'; else, ic_str = ''; end
if strcmpi(ic_mode, 'manifold'), ic_str = [ic_str ', $\mathbf{e}(0) = \mathbf{0}$']; end
if show_title
    title(ax, {sprintf('%s - %s, $\\mathbf{b}_%d \\to \\mathbf{c}_%d$, $k_s = %.2f$', osc_names{OsID}, ...
               motif_names{topology_idx}, find(bs), find(cs), ks_val), ...
           }, ...
        'Interpreter', 'latex', 'FontSize', font_title);   % Two lines: the title fits the figure width
end

% =========================================================================
% MC STEP 5b: CONVERGENCE OF THE INDEX (post-processing of the valid runs)
%   The index and its uncertainty are recomputed for every number of runs
%   n = n_start, ..., N (first n valid runs, in run order; N = n_valid).
%   (a) Deviation Index(n) - Index(N) of every non-zero coefficient,
%       coloured by mismatch type (a, b, c), inside the envelope
%       +/- 2 x max uncertainty(n) (about 95%): the curves collapse onto 0 dB as n grows.
%   (b) Largest and median uncertainty over the coefficients vs n (log-log)
%       with the reference slope n^(-1/2) of a converging Monte Carlo estimate.
% =========================================================================
active_final = S_mean > zero_rel_tol * max(S_mean);
n_start  = min(3, n_valid);                  % First n shown (n >= 2 needed for a spread; 3 for a stable one)
n_curve  = (n_start:n_valid).';
idx_n    = NaN(numel(n_curve), Ncoef);
unc_n    = NaN(numel(n_curve), Ncoef);
for m = 1:numel(n_curve)
    [idx_n(m, :), unc_n(m, :)] = indexStat(R2(1:n_curve(m), :), stat_predictor);
end
dev_n    = idx_n(:, active_final) - idx_n(end, active_final);      % deviation from the final index
umax_n   = max(unc_n(:, active_final), [], 2);
umed_n   = median(unc_n(:, active_final), 2);

fig_conv = figure('Name', 'Monte Carlo convergence', 'Position', [100 100 900 750]);   % Room for labels and legends
tlc = tiledlayout(fig_conv, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) Deviation from the final index
axc = nexttile(tlc);
hold(axc, 'on');
par_type  = mod(ceil((1:Ncoef) / numel(err_idx)) - 1, 3) + 1;      % 1 = a, 2 = b, 3 = c
par_type  = par_type(active_final);
type_col  = [0 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188];
type_name = cellfun(@(p) sprintf('$\\delta %s$', p), par_sym, 'UniformOutput', false);
h_type = gobjects(1, 3);
for tp = 1:3
    cols_tp = find(par_type == tp);
    if isempty(cols_tp), continue; end
    h = plot(axc, n_curve, dev_n(:, cols_tp), '-', 'Color', type_col(tp, :), 'LineWidth', 0.8);
    h_type(tp) = h(1);
end
h_env = plot(axc, n_curve, 2*umax_n, 'k--', 'LineWidth', 1.5);     % +/- 2 x uncertainty (about 95%)
plot(axc, n_curve, -2*umax_n, 'k--', 'LineWidth', 1.5);
yline(axc, 0, 'k-', 'LineWidth', 0.5);
grid(axc, 'on');
set(axc, 'FontSize', font_ticks, 'TickLabelInterpreter', 'latex');
xticklabels(axc, {});
ylabel(axc, 'Deviation (dB)', 'Interpreter', 'latex', 'FontSize', font_ylabel);   % Index(n) - Index(N)
y_lim = 1.1 * max(max(2 * umax_n), max(abs(dev_n(:))));          % Whole +/- 2 x uncertainty envelope visible
if y_lim > 0
    ylim(axc, [-y_lim, y_lim]);
end
ok_type = isgraphics(h_type);
legend(axc, [h_type(ok_type), h_env], [type_name(ok_type), {'$\pm 2\times$ max uncertainty'}], ...
    'Interpreter', 'latex', 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', font_ticks);
if show_title
    title(axc, sprintf('Convergence of the %s index (%d coefficients, $N_M = %d$ realizations)', ...
        lower(stat_predictor), nnz(active_final), n_valid), 'Interpreter', 'latex', 'FontSize', font_title);
end

% (b) Uncertainty vs number of runs
axd = nexttile(tlc);
loglog(axd, n_curve, umax_n, 'k-', 'LineWidth', 1.5);
hold(axd, 'on');
loglog(axd, n_curve, umed_n, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
loglog(axd, n_curve, umax_n(end) * sqrt(n_valid ./ n_curve), 'r:', 'LineWidth', 1.5);
grid(axd, 'on');
set(axd, 'FontSize', font_ticks, 'TickLabelInterpreter', 'latex');
ylabel(axd, 'Uncertainty (dB)', 'Interpreter', 'latex', 'FontSize', font_ylabel);
y_rng   = [0.8 * min(umed_n), 1.25 * max(umax_n)];
nice    = [0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100];
y_ticks = nice(nice >= y_rng(1) & nice <= y_rng(2));
if all(isfinite(y_rng)) && y_rng(2) > y_rng(1)
    set(axd, 'YLim', y_rng);
    if numel(y_ticks) >= 2
        set(axd, 'YTick', y_ticks, 'YTickLabel', compose('%g', y_ticks));
    end
end
legend(axd, {'max', 'median', '$\propto n^{-1/2}$'}, ...
    'Interpreter', 'latex', 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', font_ticks);
xlabel(tlc, 'Number of realizations $n$', 'Interpreter', 'latex', 'FontSize', font_xlabel);
x_ticks = unique([n_start, 5, 10, 20, 50, 100, 200, 500, n_valid]);
x_ticks = x_ticks(x_ticks >= n_start & x_ticks <= n_valid);
set([axc, axd], 'XScale', 'log', 'XLim', [n_start n_valid], 'XTick', x_ticks);
xticklabels(axc, {});
if ~isempty(r_conv) && r_conv < n_valid           % Where the stopping test was first met
    xline(axc, r_conv, 'k-.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    xline(axd, r_conv, 'k-.', sprintf('test met (%d runs)', r_conv), 'LineWidth', 1.0, ...
        'Interpreter', 'latex', 'FontSize', font_ticks, 'LabelVerticalAlignment', 'bottom', ...
        'HandleVisibility', 'off');
end

% =========================================================================
% MC STEP 5c: WORST RUN vs TYPICAL RUN (plot_worst_run)
%   Worst run: largest excess of its dB value over the median, over all
%   non-zero coefficients. Typical run: closest to the median in that
%   coefficient. Both are integrated again (same seed, same solver) and the
%   stored RMS is recomputed as a check. Burst diagnostic: time the
%   coefficient stays above 10 x the median RMS; a burst lasting many steps
%   (many dt) is resolved dynamics, one lasting a few steps is numerical.
% =========================================================================
worst_info = [];
fig_w = [];
if plot_worst_run && n_valid > 1
    D_all  = 20*log10(RMS_runs);                          % r_done x 54, NaN for blown runs
    med_dB = median(D_all, 1, 'omitnan');
    excess = max(D_all, [], 1, 'omitnan') - med_dB;
    excess(zero_coef) = -Inf;
    [exc_w, k_w] = max(excess);                           % coefficient
    [~, r_w]     = max(D_all(:, k_w));                    % worst run
    [~, r_m]     = min(abs(D_all(:, k_w) - med_dB(k_w))); % typical run
    j_w = ceil(k_w / numel(err_idx));                     % mismatch q_j
    i_w = k_w - (j_w - 1) * numel(err_idx);               % error component
    fprintf('\nWorst run: run %d (seed %d), %s = %.2f dB, %.2f dB above the median (%.2f dB; typical run %d).\n', ...
        r_w, base_seed + r_w, S_labels{k_w}, D_all(r_w, k_w), exc_w, med_dB(k_w), r_m);

    rows_w = [1:n, err_idx, sens_rows(k_w)];              % x1, errors, s_k
    runs_w = [r_w, r_m];
    T_w = cell(1, 2);
    Y_w = cell(1, 2);
    for q = 1:2
        rr = runs_w(q);
        rng(base_seed + rr);                              % Same initial condition as in the loop
        X0 = reshape(rand(n * num_oscillators, 1), n, num_oscillators);
        E0 = [X0(:,1) - X0(:,2); X0(:,1) - X0(:,3); X0(:,2) - X0(:,3)];
        if strcmpi(ic_mode, 'manifold')
            E0 = zeros(size(E0));                  % Start on the synchronization manifold
        end
        if use_pool
            x1_0 = X1_pool(rr, :).';
        else
            x1_0 = X0(:,1);
        end
        Z0 = zeros(Nz, 1);
        Z0(1:Nx) = [x1_0; E0];
        [T_w{q}, Y_w{q}] = integrateRun(Z0, solver, f_rhs, f_ode, ode_opts, dt, N_steps, rows_w);
    end

    % Check: RMS recomputed on the same samples as in the loop
    idx_s = (k_trans + k_sample : k_sample : N_steps) + 1;
    sb_w  = abs(Y_w{1}(end, :)) * sens_scal(k_w);         % |sbar_k(t)|, worst run
    sb_m  = abs(Y_w{2}(end, :)) * sens_scal(k_w);         % |sbar_k(t)|, typical run
    fprintf('  Recomputed RMS: worst %.2f dB (stored %.2f), typical %.2f dB (stored %.2f)\n', ...
        20*log10(sqrt(mean(sb_w(idx_s).^2))), D_all(r_w, k_w), ...
        20*log10(sqrt(mean(sb_m(idx_s).^2))), D_all(r_m, k_w));

    % Burst diagnostic in the window t >= t_transient
    thr    = 10 * 10^(med_dB(k_w) / 20);                  % 10 x median RMS
    in_win = T_w{1} >= t_transient;
    above  = in_win & (sb_w > thr);
    d      = diff([0, above, 0]);
    b_beg  = find(d == 1);
    b_end  = find(d == -1) - 1;
    dT     = mean(diff(T_w{1}));
    if isempty(b_beg)
        fprintf('  Never above 10 x median RMS in the window.\n');
        burst_len = 0;
    else
        burst_len = (b_end - b_beg + 1) * dT;             % duration of each burst
        fprintf(['  Above 10 x median RMS: %d burst(s), %.1f time units in total (%.1f%% of the window);\n' ...
                 '  longest burst %.2f time units = %d steps of %g\n'], ...
            numel(b_beg), sum(burst_len), 100 * sum(burst_len) / (T_end - t_transient), ...
            max(burst_len), round(max(burst_len) / dt), dt);
    end
    worst_info = struct('run', r_w, 'typical_run', r_m, 'coef', S_labels{k_w}, ...
        'dB', D_all(r_w, k_w), 'median_dB', med_dB(k_w), 'burst_len', burst_len);

    % Figure: |sbar_k|, ||e|| and node 1 for both runs
    fig_w = figure('Name', 'Worst run vs typical run');
    tlw = tiledlayout(fig_w, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    axw1 = nexttile(tlw);
    semilogy(axw1, T_w{2}, sb_m, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);
    hold(axw1, 'on');
    semilogy(axw1, T_w{1}, sb_w, 'r', 'LineWidth', 1.0);
    yline(axw1, 10^(med_dB(k_w) / 20), 'k--', 'LineWidth', 1.0);
    xline(axw1, t_transient, 'k:', 'LineWidth', 1.0);
    grid(axw1, 'on');
    set(axw1, 'FontSize', font_ticks, 'TickLabelInterpreter', 'latex');
    xticklabels(axw1, {});
    ylabel(axw1, sprintf('$|\\bar{s}|$: $\\partial %s / \\partial %s$', comp_tex{i_w}, par_tex{j_w}), ...
        'Interpreter', 'latex', 'FontSize', font_ylabel);
    legend(axw1, {sprintf('Typical run %d', r_m), sprintf('Worst run %d', r_w), 'Median RMS'}, ...
        'Interpreter', 'latex', 'Location', 'best');

    axw2 = nexttile(tlw);
    semilogy(axw2, T_w{2}, vecnorm(Y_w{2}(n+1:n+numel(err_idx), :), 2, 1), 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);
    hold(axw2, 'on');
    semilogy(axw2, T_w{1}, vecnorm(Y_w{1}(n+1:n+numel(err_idx), :), 2, 1), 'r', 'LineWidth', 1.0);
    xline(axw2, t_transient, 'k:', 'LineWidth', 1.0);
    grid(axw2, 'on');
    set(axw2, 'FontSize', font_ticks, 'TickLabelInterpreter', 'latex');
    xticklabels(axw2, {});
    ylabel(axw2, '$\|\mathbf{e}\|$', 'Interpreter', 'latex', 'FontSize', font_ylabel);

    axw3 = nexttile(tlw);
    plot(axw3, T_w{1}, Y_w{1}(1, :), 'b', 'LineWidth', 1.0);
    hold(axw3, 'on');
    plot(axw3, T_w{1}, Y_w{1}(3, :), 'Color', [0 0.6 0], 'LineWidth', 1.0);
    grid(axw3, 'on');
    set(axw3, 'FontSize', font_ticks, 'TickLabelInterpreter', 'latex');
    ylabel(axw3, sprintf('Node 1 (run %d)', r_w), 'Interpreter', 'latex', 'FontSize', font_ylabel);
    legend(axw3, {'$x_1$', '$z_1$'}, 'Interpreter', 'latex', 'Location', 'best');

    xlabel(tlw, '$t$', 'Interpreter', 'latex', 'FontSize', font_xlabel);
    title(tlw, sprintf('%s - %s: worst run %d vs typical run %d (%s)', osc_names{OsID}, ...
        motif_names{topology_idx}, r_w, r_m, solver), 'Interpreter', 'latex', 'FontSize', font_title);
    linkaxes([axw1, axw2, axw3], 'x');
    xlim(axw1, [0 T_end]);
end

% Hide the axes toolbars so they are not captured by exportgraphics
all_ax = findall([fig, fig_conv, fig_w], 'Type', 'axes');
for kk = 1:numel(all_ax)
    all_ax(kk).Toolbar.Visible = 'off';
end

% =========================================================================
% MC STEP 6: SAVE (mat/ and eps/)
% =========================================================================
file_tag = sprintf('%s_%s_B%d_ks%.2f_%s', FMtype, topology_list{topology_idx}, CasoBs, ks_val, stat_predictor);
if use_pool
    file_tag = [file_tag '_pool'];                    % Separate files to compare both options
end
if strcmpi(solver, 'ode45')
    file_tag = [file_tag '_ode45'];                   % RK4 file names unchanged
end
if strcmpi(ic_mode, 'manifold')
    file_tag = [file_tag '_manifold'];                % Previous (random e0) files are kept
end
if ~exist('eps', 'dir'), mkdir('eps'); end
if ~exist('mat', 'dir'), mkdir('mat'); end
exportgraphics(fig, fullfile('eps', sprintf('Motif_SensitivityIndex_dB_%s.eps', file_tag)), ...
    'ContentType', 'vector');
exportgraphics(fig_conv, fullfile('eps', sprintf('Motif_SensitivityIndex_Convergence_%s.eps', file_tag)), ...
    'ContentType', 'vector');
save(fullfile('mat', sprintf('SensitivityIndex_%s.mat', file_tag)), ...
    'S_table', 'RMS_runs', 'S_dB_map', 'Index_dB', 'Unc_dB', 'growth_run', 'Parval', 'OscPar', ...
    'p_scale', 'xref', 'T_end', 'dt', 'dt_sample', 't_transient', 'r_done', 'n_valid', 'blown_run', ...
    'converged', 'r_conv', 'stop_on_convergence', 'n_curve', 'idx_n', 'unc_n', 'hist_r', 'hist_idx', 'delta_hist', 'semdB_hist', 'stat_predictor', ...
    'err_bound', 'conv_tol_dB', 'conv_sem_dB', 'conv_window', 'zero_rel_tol', ...
    'base_seed', 'OsID', 'topology_idx', 'CasoBs', 'ks_val', ...
    'use_pool', 'X1_pool', 'pool_T_transient', 'pool_spacing', 'ic_mode', ...
    'solver', 'ode_RelTol', 'ode_AbsTol', 'worst_info', 'V_table', 'V_table_dB', 'V_grid', 'V_grid_unc');

% LaTeX table of the sensitivity vector index (tex/)
if ~exist('tex', 'dir'), mkdir('tex'); end
tex_file = fullfile('tex', sprintf('SensitivityVector_%s.tex', file_tag));
fid = fopen(tex_file, 'w');
fprintf(fid, '%% RMS sensitivity vector index: %s - %s, b_%d -> c_%d, ks = %.2f, %d runs, %s\n', ...
    FMtype, topology_list{topology_idx}, find(bs), find(cs), ks_val, n_valid, stat_predictor);
fprintf(fid, '%% All values in dB; the norm is computed on linear values, then 20 log10\n');
fprintf(fid, '\\begin{tabular}{cccccccc}\n\\toprule\n');
fprintf(fid, ['\\textbf{Sensitivity vector} & \\textbf{Parameter} & \\textbf{Error vector} & ' ...
              '$e_x$ (dB) & $e_y$ (dB) & $e_z$ (dB) & $\\|s_{q_j}(\\mathbf{e})\\|$ (dB) & $\\pm$ \\textbf{Uncert.} (dB) \\\\\n']);
fprintf(fid, '\\midrule\n');
for row = 1:nV
    par_t = sprintf('%s_%s', V_param{row}(1), V_param{row}(2));        % a_2
    err_t = sprintf('\\mathbf{e}_{%s}', V_err{row}(2:3));              % \mathbf{e}_{12}
    if V_norm(row) == 0
        unc_t = '--';                                                  % Zero vector: no uncertainty
    else
        unc_t = dbStr(V_unc_dB(row), true);
    end
    fprintf(fid, '$s_{%s}(%s)$ & $%s$ & $%s$ & %s & %s & %s & %s & %s \\\\\n', ...
        par_t, err_t, par_t, err_t, dbStr(V_comp_dB(row, 1), true), dbStr(V_comp_dB(row, 2), true), ...
        dbStr(V_comp_dB(row, 3), true), dbStr(V_norm_dB(row), true), unc_t);
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);
fprintf('Saved: %s\n', tex_file);

% Alternative LaTeX table, heat-map layout (kept for comparison)
tex_file_grid = fullfile('tex', sprintf('SensitivityVectorGrid_%s.tex', file_tag));
fid = fopen(tex_file_grid, 'w');
fprintf(fid, '%% RMS sensitivity vector index: %s - %s, b_%d -> c_%d, ks = %.2f, %d runs, %s\n', ...
    FMtype, topology_list{topology_idx}, find(bs), find(cs), ks_val, n_valid, stat_predictor);
fprintf(fid, '%% Same layout as the heat map; components and norms in linear units\n');
fprintf(fid, '\\begin{tabular}{ll%s}\n\\toprule\n', repmat('c', 1, Nq));
fprintf(fid, '\\textbf{Error vector} & \\textbf{Comp.}');
for j = 1:Nq
    fprintf(fid, ' & $%s$', par_tex{j});
end
fprintf(fid, ' \\\\\n\\midrule\n');
comp_tex_v = {'$x$', '$y$', '$z$', '$\|s_{q_j}(\mathbf{e})\|$'};   % printed with %s: single backslashes
for p = 1:num_pairs
    for c = 1:n + 1
        rg = (p - 1) * (n + 1) + c;
        if c == 1
            fprintf(fid, '$\\mathbf{e}_{%s}$', pair_lbl{p}(2:3));
        end
        fprintf(fid, ' & %s', comp_tex_v{c});
        for j = 1:Nq
            if c <= n
                fprintf(fid, ' & %.4g', V_grid(rg, j));
            elseif V_grid(rg, j) == 0
                fprintf(fid, ' & 0');
            else
                fprintf(fid, ' & $%.4g \\pm %.2g$', V_grid(rg, j), V_grid_unc(rg, j));
            end
        end
        fprintf(fid, ' \\\\\n');
    end
    if p < num_pairs, fprintf(fid, '\\midrule\n'); end
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);
fprintf('Saved: %s\n', tex_file_grid);
if ~isempty(fig_w)
    exportgraphics(fig_w, fullfile('eps', sprintf('Motif_SensitivityIndex_WorstRun_%s.eps', file_tag)), ...
        'ContentType', 'vector');
end
fprintf('Saved: eps/Motif_SensitivityIndex_dB_%s.eps and mat/SensitivityIndex_%s.mat\n', file_tag, file_tag);

% =========================================================================
% LOCAL FUNCTION: INDEX AND ITS UNCERTAINTY (dB) FROM PER-RUN RMS VALUES
%   R: runs x coefficients (linear RMS)
% =========================================================================
function [idx_dB, unc_dB] = indexStat(R, predictor)
    nr = size(R, 1);
    switch lower(predictor)
        case 'sem'
            m      = mean(R, 1);
            idx_dB = 20*log10(m);
            unc_dB = 20*log10(1 + std(R, 0, 1) ./ sqrt(nr) ./ m);
        case 'median'
            D      = 20*log10(R);
            idx_dB = median(D, 1);
            unc_dB = 1.2533 * std(D, 0, 1) ./ sqrt(nr);   % standard error of the median
        otherwise
            error('stat_predictor must be ''SEM'' or ''Median''.');
    end
    if nr < 2
        unc_dB(:) = Inf;
    end
end

% =========================================================================
% LOCAL FUNCTION: ONE RUN WITH THE FULL TIME SERIES (worst-run check)
%   Output every dt on [0, T_end]; only the requested rows are kept.
% =========================================================================
function [t, Y] = integrateRun(Z0, solver, f_rhs, f_ode, ode_opts, dt, N_steps, rows)
    if strcmpi(solver, 'ode45')
        [t, Z] = ode45(f_ode, (0:N_steps) * dt, Z0, ode_opts);
        t = t.';
        Y = Z(:, rows).';
    else
        t = (0:N_steps) * dt;
        Y = zeros(numel(rows), N_steps + 1);
        Y(:, 1) = Z0(rows);
        Z = Z0;
        for k = 1:N_steps
            k1 = f_rhs(Z);
            k2 = f_rhs(Z + 0.5*dt*k1);
            k3 = f_rhs(Z + 0.5*dt*k2);
            k4 = f_rhs(Z + dt*k3);
            Z  = Z + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
            Y(:, k+1) = Z(rows);
        end
    end
end

% =========================================================================
% LOCAL FUNCTION: BLOW-UP EVENT FOR ode45
%   Stops the run when any |e_ij| exceeds err_bound or the state is NaN/Inf.
% =========================================================================
function [value, isterminal, direction] = blowupEvent(Z, err_idx, err_bound)
    if any(~isfinite(Z))
        value = -1;
    else
        value = err_bound - max(abs(Z(err_idx)));
    end
    isterminal = 1;
    direction  = 0;
end

% =========================================================================
% LOCAL FUNCTION: dB VALUE AS TEXT (2 decimals, -Inf for exact zeros,
%   no '-0.00'); as_tex = true writes LaTeX ($-\infty$)
% =========================================================================
function str = dbStr(x, as_tex)
    if isinf(x) && x < 0
        if as_tex, str = '$-\infty$'; else, str = '-Inf'; end
    else
        str = sprintf('%.2f', x);
        if strcmp(str, '-0.00'), str = '0.00'; end
    end
end

% =========================================================================
% LOCAL FUNCTION: STRING CHOICE (inline if/else for the table print)
% =========================================================================
function str = ternary_str(cond, str_true, str_false)
    if cond
        str = str_true;
    else
        str = str_false;
    end
end

% =========================================================================
% LOCAL FUNCTION: COLUMN-BY-COLUMN EVALUATION (fallback if not vectorized)
% =========================================================================
function F = evalColumns(f_one, Z)
    F = zeros(size(Z));
    for col = 1:size(Z, 2)
        F(:, col) = f_one(Z(:, col));
    end
end

% =========================================================================
% LOCAL FUNCTION: PERCENTILE PER COLUMN (no toolbox needed, NaN ignored)
% =========================================================================
function v = colPrctile(X, pct)
    v = NaN(1, size(X, 2));
    for col = 1:size(X, 2)
        x = sort(X(~isnan(X(:, col)), col));
        m = numel(x);
        if m == 0, continue; end
        pos = 1 + (m - 1) * pct / 100;          % Linear interpolation
        lo  = floor(pos);
        hi  = min(lo + 1, m);
        v(col) = x(lo) + (pos - lo) * (x(hi) - x(lo));
    end
end
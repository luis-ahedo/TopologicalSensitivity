% =========================================================================
% SCRIPT: MotifLLE_vs_ks_Case.m
% DESCRIPTION: Sweep of the coupling strength ks in I_ks = [0, kf] for one
%              topology and coupling scheme. For each ks the largest
%              transverse Lyapunov exponent of the error system (linearized
%              on e = 0) is computed with MaxTransLyaExp; the result is
%              plotted as LLE(ks) vs ks. The symbolic part generates
%              odeMotifVariational.m.
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
n = length(X);          % State dimension per node (3)
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
[bs, cs] = SISObc(CasoBs, n);         % Inner coupling vectors
% =========================================================================
% STEP 3: AUTOMATED STATE & VECTOR FIELD GENERATION
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
% ERROR SYSTEM:  de_ij/dt = [F(X_i) - F(X_i - e_ij)] + u_i - u_j
%   u_k = ks * sum_{s -> k} B*C*(X_s - X_k)  is the coupling input of node k,
%   and every difference X_s - X_k is written with its own error symbol
%   (+e_sk if s < k, -e_ks if s > k). Each equation keeps e12, e13, e23;
%   the nonlinear part depends on e_ij and the reference node X_i only.
% =========================================================================
fprintf('Building the error system de12/dt, de13/dt, de23/dt...\n');

% Pair index lookup: P(i,j) = position of e_ij in E_array
P = zeros(num_oscillators);
for p = 1:num_pairs
    P(pair_ij(p,1), pair_ij(p,2)) = p;
end
errDiff = @(s, k) (s < k) * E_array{P(min(s,k), max(s,k))} ...
                - (s > k) * E_array{P(min(s,k), max(s,k))};   % X_s - X_k

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
    FXi_minus = subs(F, X, Xi - E_array{p});               % F(X_i - e_ij)
    dEij = (FXj(idx{i}) - FXi_minus) + U{i} - U{j};
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
% VARIATIONAL EQUATIONS ON THE SYNCHRONIZATION MANIFOLD (e = 0)
%   d(delta e)/dt = J_E(x1) * delta e,   J_E = d(F_matrix_E)/de at e = 0
%   dx1/dt        = F(x1)                (all coupling vanishes on e = 0)
% =========================================================================
fprintf('\nBuilding the variational equations (Jacobian at e = 0)...\n');

% Plain symbols (jacobian w.r.t. x1(t), e12x(t), ... is avoided)
xs = sym('xs', [n 1], 'real');             % node 1 states
dim_V = n * (num_oscillators - 1);         % independent errors [e12; e13]
es = sym('es', [dim_Y 1], 'real');         % errors [e12; e13; e23]
dY = sym('dY', [dim_V 1], 'real');         % perturbations [de12; de13]

% Error system in plain symbols, with X2 = X1 - e12 (used in de23/dt)
F_E_plain = subs(F_matrix_E, [Xj(idx{1}); Xj(idx{2}); vars_E], ...
                             [xs; xs - es(1:n); es]);

% Jacobian w.r.t. the errors, evaluated on the synchronization manifold
% Independent system: de12/dt and de13/dt with e23 = e13 - e12 in both.
% de23/dt is redundant; keeping it adds a fictitious-node mode whose
% exponent is not physical (e.g. the uncoupled one in the auxiliary motif).
F_E_red = subs(F_E_plain(1:dim_V), es(2*n+1:3*n), es(n+1:2*n) - es(1:n));
J_E  = jacobian(F_E_red, es(1:dim_V));
J_E0 = simplify(subs(J_E, es(1:dim_V), zeros(dim_V, 1)));   % depends on xs and Par only

% Augmented system [x1; delta e]
F_X1  = subs(FXj(idx{1}), Xj(idx{1}), xs);         % dx1/dt = F(x1)
F_var = [F_X1; J_E0 * dY];

fprintf('\nVariational equations (%s, %s, CasoBs = %d):\n', FMtype, edges_input, CasoBs);
for p = 1:dim_V/n
    fprintf('d(delta %s)/dt =\n', pair_names{p});
    disp(F_var(n + (p-1)*n + (1:n)));
end

% =========================================================================
% SINGLE ODE FUNCTION FOR [x1; delta e]  (odeMotifVariational.m)
%   Same file and signature used by MaxTransLyaExp(Parval, dim_Y, n, opts):
%   dZ = odeMotifVariational(t, Z, Parval), Z = [x1; dY] (n + dim_V states).
%   Generated code indexes Z(k,:), so Z can hold several columns
%   [x1; q_i] and all of them are evaluated in one call.
% =========================================================================
Zf     = str2sym(compose('Z%d(t)', (1:n+dim_V).'));   % state functions Z1(t)...
F_varZ = subs(F_var, [xs; dY], Zf);
filename_ode = 'odeMotifVariational.m';
fprintf('Generating MATLAB function file: %s...\n', filename_ode);
odeFunction(F_varZ, Zf, Par, 'File', filename_ode, 'Optimize', true);
clear odeMotifVariational;
rehash;
% =========================================================================
% KS SWEEP CONFIGURATION
% =========================================================================
kf     = 25;                    % I_ks = [0, kf]
N_ks   = 101;                   % Number of ks points
ks_vec = linspace(0, kf, N_ks);

rng(42);
initial_conditions = rand(total_states, 1);         % Same seed as SimulateMotif.m

opts = struct();
opts.tf          = 2000;                            % Maximum averaging time
opts.dt          = 0.005;                           % RK4 fixed step
opts.T_transient = 100;                             % Transient of node 1
opts.num_LE      = 1;                               % Largest exponent only
opts.x0          = initial_conditions(idx{1});      % Node 1 initial condition
opts.seed        = 42;                              % Same perturbation for every ks
opts.early_exit  = true;                            % Stop once LLE has converged

use_parallel = false;           % true: parfor (Parallel Computing Toolbox)

font_title  = 13;
font_xlabel = 11;
font_ylabel = 11;
font_ticks  = 10;

% =========================================================================
% SWEEP
% =========================================================================
LLE_vec = NaN(1, N_ks);
fprintf('\nSweeping ks in [0, %g] (%d points)...\n', kf, N_ks);
t_start = tic;
if use_parallel
    parfor k_idx = 1:N_ks
        LLE_vec(k_idx) = MaxTransLyaExp([OscPar, ks_vec(k_idx)], dim_V, n, opts);
    end
else
    for k_idx = 1:N_ks
        LLE_vec(k_idx) = MaxTransLyaExp([OscPar, ks_vec(k_idx)], dim_V, n, opts);
        fprintf('  ks = %6.3f  ->  LLE = %+.5f\n', ks_vec(k_idx), LLE_vec(k_idx));
    end
end
fprintf('Sweep finished in %.1f min\n', toc(t_start)/60);

% =========================================================================
% PLOT: LLE(ks) vs ks
% =========================================================================
fig_lle = figure('Name', sprintf('LLE vs ks (%s - %s, CasoBs = %d)', FMtype, edges_input, CasoBs));
ax = axes(fig_lle);
plot(ax, ks_vec, LLE_vec, 'b-', 'LineWidth', 1.5);
hold(ax, 'on');
yline(ax, 0, 'k--', 'LineWidth', 1.0);
grid(ax, 'on');
xlim(ax, [0 kf]);
set(ax, 'FontSize', font_ticks);
xlabel(ax, '$k_s$', 'Interpreter', 'latex', 'FontSize', font_xlabel);
ylabel(ax, '$\lambda^\perp_{\max}$', 'Interpreter', 'latex', 'FontSize', font_ylabel);
[bs_t, cs_t] = SISObc(CasoBs, n);
title(ax, sprintf('Largest Transverse Lyapunov Exponent ($\\mathbf{b}_%d \\to \\mathbf{c}_%d$)', find(bs_t), find(cs_t)), ...
    'Interpreter', 'latex', 'FontSize', font_title);

% =========================================================================
% EPS EXPORT
% =========================================================================
folder_name = 'eps';
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end
file_tag = sprintf('%s_%s_B%d_kf%g', FMtype, edges_input, CasoBs, kf);
filename_eps = fullfile(folder_name, sprintf('Motif_LLE_vs_ks_%s.eps', file_tag));
exportgraphics(fig_lle, filename_eps, 'ContentType', 'vector');
fprintf('Saved: %s\n', filename_eps);

% =========================================================================
% SCRIPT: MotifErrorSystem.m
% DESCRIPTION: Symbolic error system of a 3-node chaotic network motif,
%              e_ij = X_i - X_j, as three equations de12/dt, de13/dt,
%              de23/dt. Coupling terms are written with their own error
%              symbols; nonlinear terms as F(X_i) - F(X_i - e_ij).
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
ks_val = 2.90;                        % Scalar coupling strength
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
% SUBSTITUTION OF e23 IN de12/dt AND de13/dt:  e23 = e13 - e12
%   (rows 1:6, every topology; de12/dt contains e23 when node 2 is driven
%   by node 3, e.g. relay, mixed master, mixed slave and the rings)
% =========================================================================
F_matrix_E(1:6) = expand(subs(F_matrix_E(1:6), E_array{3}, E_array{2} - E_array{1}));
for r = 1:6
    F_matrix_E(r) = collect(F_matrix_E(r), Par);
end
fprintf('\nde12/dt after substituting e23 = e13 - e12:\n');
disp(F_matrix_E(1:3));
fprintf('\nde13/dt after substituting e23 = e13 - e12:\n');
disp(F_matrix_E(4:6));

% =========================================================================
% AUGMENTED CLOSED SYSTEM [X1; e12; e13] AND ODE FUNCTION
%   Node 1 supplies the reference trajectory. de23/dt is redundant and is
%   not integrated: e23 = e13 - e12 is rebuilt after the integration.
% =========================================================================
X1v      = Xj(idx{1});
F_X1     = FXj(idx{1}) + U{1};                              % dX1/dt, coupling in errors
F_aug    = [F_X1; F_matrix_E(1:6)];
F_aug    = subs(F_aug, Xj(idx{2}), X1v - E_array{1});       % X2 = X1 - e12
vars_aug = [X1v; vars_E(1:6)];

odeFunction(F_aug, vars_aug, Par, 'File', 'odeMotifError.m', 'Optimize', true);
clear odeMotifError;
rehash;

% =========================================================================
% NUMERICAL INTEGRATION
% =========================================================================
rng(42);
initial_conditions = rand(total_states, 1);                 % Same node ICs as SimulateMotif.m
X0 = reshape(initial_conditions, n, num_oscillators);       % Column k = node k
Y0 = [X0(:,1); X0(:,1) - X0(:,2); X0(:,1) - X0(:,3)];

tspan  = 0:0.01:100;
Parval = [OscPar(1), OscPar(2), OscPar(3), ks_val];
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

fprintf('Integrating the error system (ks = %.3f)...\n', ks_val);
[t_out, Y_out] = ode45(@(t, Y) odeMotifError(t, Y, Parval), tspan, Y0, options);

E_out = [Y_out(:, n+1:3*n), Y_out(:, 2*n+1:3*n) - Y_out(:, n+1:2*n)];   % [e12 e13 e23], 9 columns

% =========================================================================
% PLOT: ERROR COMPONENTS OF e12, e13, e23
% =========================================================================
font_title  = 20;   % Global title (tiledlayout)
font_xlabel = 20;
font_ylabel = 20;
font_legend = 15;
font_ticks  = 15;

state_letters = {'x', 'y', 'z'};
fig_err = figure('Name', 'Error System');
tl = tiledlayout(num_pairs + 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Node 1 dynamics (reference trajectory)
ax = nexttile(tl);
plot(ax, t_out, Y_out(:, 1:n), 'LineWidth', 1.2);
grid(ax, 'on');
set(ax, 'FontSize', font_ticks);
ylabel(ax, '$\mathbf{x}_1$', 'Interpreter', 'latex', 'FontSize', font_ylabel);
legend(ax, cellfun(@(s) sprintf('$%s_1$', s), state_letters, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'FontSize', font_legend, 'Location', 'northeast', 'Orientation', 'horizontal');

for p = 1:num_pairs
    ax = nexttile(tl);
    plot(ax, t_out, E_out(:, (p-1)*n + (1:n)), 'LineWidth', 1.2);
    grid(ax, 'on');
    set(ax, 'FontSize', font_ticks);
    ylabel(ax, sprintf('$\\mathbf{e}_{%s}$', pair_names{p}(2:3)), 'Interpreter', 'latex', 'FontSize', font_ylabel);
    legend(ax, cellfun(@(s) sprintf('$e_{%s%s}$', pair_names{p}(2:3), s), state_letters, 'UniformOutput', false), ...
        'Interpreter', 'latex', 'FontSize', font_legend, 'Location', 'northeast', 'Orientation', 'horizontal');
end
xlabel(tl, 'Time $t$', 'Interpreter', 'latex', 'FontSize', font_xlabel);
title(tl, sprintf('Error System: %s ($k_s = %.2f$)', motif_names{topology_idx}, ks_val), ...
    'Interpreter', 'latex', 'FontSize', font_title);

% =========================================================================
% EPS EXPORT
% =========================================================================
folder_name = 'eps';
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end
file_tag = sprintf('%s_%s_B%d_ks%.2f', FMtype, edges_input, CasoBs, ks_val);
filename_eps = fullfile(folder_name, sprintf('Motif_ErrorSystem_%s.eps', file_tag));
exportgraphics(fig_err, filename_eps, 'ContentType', 'vector');
fprintf('Saved: %s\n', filename_eps);

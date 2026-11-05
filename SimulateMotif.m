% =========================================================================
% SCRIPT: SimulateMotif.m
% DESCRIPTION: Full simulation pipeline for 3-node chaotic network motifs.
%              Includes symbolic setup, global coupling matrix assembly,
%              standalone ODE generation, numerical simulation with blow-up
%              protection, tiled trajectory/error visualization, manifold
%              projections, 3D phase space, and automated vector EPS export.
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
total_states = n * num_oscillators;
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

F_coupled = FXj + Mc * Xj;

% =========================================================================
% STEP 5: AUTOMATED ODE FUNCTION GENERATION
% =========================================================================
filename = 'odeMotifCoupled.m';
fprintf('Generating MATLAB function file: %s...\n', filename);
odeFunction(F_coupled, Xj, Par, 'File', filename, 'Optimize', true);
clear odeMotifCoupled;   % avoid reusing a cached old version of the file
rehash;

% =========================================================================
% STEP 6: NUMERICAL SIMULATION & DATA EXTRACTION
% =========================================================================
font_title  = 13;     
font_xlabel = 11;     
font_ylabel = 11;     
font_zlabel = 11;
font_legend = 10;     
font_ticks  = 10;     

Parval = [OscPar(1), OscPar(2), OscPar(3), ks_val];

rng(42);
initial_conditions = rand(total_states, 1);  
tspan = [0 100];                             

% Blow-up event options
max_state_bound = 1e4;
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, ...
                 'Events', @(t, Y) blowupEvent(t, Y, max_state_bound));

fprintf('Solving the coupled network (ks = %.3f)...\n', ks_val);
[t_out, Y_out, te] = ode45(@(t, Y) odeMotifCoupled(t, Y, Parval), tspan, initial_conditions, options);
if ~isempty(te)
    warning('Blow-up detected at t = %.2f: plots stop there.', te(1));
end

x_states = zeros(length(t_out), num_oscillators);
y_states = zeros(length(t_out), num_oscillators);
z_states = zeros(length(t_out), num_oscillators);

for k = 1:num_oscillators
    x_states(:, k) = Y_out(:, (k-1)*n + 1);
    y_states(:, k) = Y_out(:, (k-1)*n + 2);
    z_states(:, k) = Y_out(:, (k-1)*n + 3);
end

pairs_matrix = nchoosek(1:num_oscillators, 2);
num_pairs = size(pairs_matrix, 1);
e_norm = zeros(length(t_out), num_pairs);
sync_err_labels = cell(num_pairs, 1);

for p_idx = 1:num_pairs
    idx_i = pairs_matrix(p_idx, 1);
    idx_j = pairs_matrix(p_idx, 2);
    ex = x_states(:, idx_i) - x_states(:, idx_j);
    ey = y_states(:, idx_i) - y_states(:, idx_j);
    ez = z_states(:, idx_i) - z_states(:, idx_j);
    e_norm(:, p_idx) = sqrt(ex.^2 + ey.^2 + ez.^2);
    sync_err_labels{p_idx} = sprintf('$e_{%d%d}(t)$', idx_i, idx_j);
end

% =========================================================================
% STEP 7: TILED LAYOUT (TRAJECTORIES & SYNCHRONIZATION ERRORS)
% =========================================================================
fig_combined = figure('Name', 'State Dynamics & Synchronization Error');
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Tile 1: State Trajectories ---
ax1 = nexttile(tl);
hold(ax1, 'on');
plot(ax1, t_out, x_states, 'LineWidth', 1.2);
plot(ax1, t_out, y_states, 'LineWidth', 1.2);
plot(ax1, t_out, z_states, 'LineWidth', 1.2);
grid(ax1, 'on');

% Dynamically build state variable labels
state_traj_labels = cell(3 * num_oscillators, 1);
lbl_cnt = 1;
for k = 1:num_oscillators
    state_traj_labels{lbl_cnt} = sprintf('$x_%d(t)$', k);
    lbl_cnt = lbl_cnt + 1;
end
for k = 1:num_oscillators
    state_traj_labels{lbl_cnt} = sprintf('$y_%d(t)$', k);
    lbl_cnt = lbl_cnt + 1;
end
for k = 1:num_oscillators
    state_traj_labels{lbl_cnt} = sprintf('$z_%d(t)$', k);
    lbl_cnt = lbl_cnt + 1;
end

set(ax1, 'FontSize', font_ticks);
xticklabels(ax1, {});  % Hide top x-ticks to share x-axis with Tile 2
ylabel(ax1, '$\mathbf{x}_k(t)$', 'Interpreter', 'latex', 'FontSize', font_ylabel);
title(ax1, sprintf('Oscillator Trajectories ($N = %d$, $k_s = %.2f$)', num_oscillators, ks_val), ...
    'Interpreter', 'latex', 'FontSize', font_title);

% Horizontal legend placed south-outside Tile 1
legend(ax1, state_traj_labels, 'Interpreter', 'latex', 'FontSize', font_legend, ...
    'Orientation', 'horizontal', 'Location', 'southoutside');

% --- Tile 2: Synchronization Errors ---
ax2 = nexttile(tl);
plot(ax2, t_out, e_norm, 'LineWidth', 1.5);
grid(ax2, 'on');

set(ax2, 'FontSize', font_ticks);
xlabel(ax2, 'Time $t$', 'Interpreter', 'latex', 'FontSize', font_xlabel);
ylabel(ax2, '$\|e_{ij}(t)\|$', 'Interpreter', 'latex', 'FontSize', font_ylabel);
title(ax2, 'Pairwise Synchronization Errors', 'Interpreter', 'latex', 'FontSize', font_title);
legend(ax2, sync_err_labels, 'Interpreter', 'latex', 'FontSize', font_legend, 'Location', 'best');

% Link X-Axes across both subplots
linkaxes([ax1, ax2], 'x');

% =========================================================================
% STEP 8: SYNCHRONIZATION MANIFOLD PROJECTIONS (1:1 PLANES)
% =========================================================================
fprintf('Rendering pairwise synchronization manifold projections...\n');

fig_planes = figure('Name', 'Synchronization Manifold Projections');
planes_layout = tiledlayout(3, num_pairs, 'TileSpacing', 'compact', 'Padding', 'compact');

state_mats    = {x_states, y_states, z_states};
state_letters = {'x', 'y', 'z'};
plane_colors  = {'b', 'g', 'm'};

for row = 1:3
    S = state_mats{row};
    for p_idx = 1:num_pairs
        idx_i = pairs_matrix(p_idx, 1);
        idx_j = pairs_matrix(p_idx, 2);
        
        ax_plane = nexttile(planes_layout);
        plot(ax_plane, S(:, idx_i), S(:, idx_j), plane_colors{row}, 'LineWidth', 1.2);
        hold(ax_plane, 'on');
        
        % 1:1 Synchronization Diagonal Line
        all_vals = [S(:, idx_i); S(:, idx_j)];
        min_v = min(all_vals); max_v = max(all_vals);
        plot(ax_plane, [min_v, max_v], [min_v, max_v], 'k--', 'LineWidth', 1.2);
        grid(ax_plane, 'on');
        
        set(ax_plane, 'FontSize', font_ticks);
        xlabel(ax_plane, sprintf('$%s_%d(t)$', state_letters{row}, idx_i), 'Interpreter', 'latex', 'FontSize', font_xlabel);
        ylabel(ax_plane, sprintf('$%s_%d(t)$', state_letters{row}, idx_j), 'Interpreter', 'latex', 'FontSize', font_ylabel);
    end
end

title(planes_layout, sprintf('Synchronization Manifold Projections ($k_s = %.2f$)', ks_val), ...
    'Interpreter', 'latex', 'FontSize', font_title);

% =========================================================================
% STEP 9: 3D PHASE-SPACE TRAJECTORIES
% =========================================================================
fprintf('Rendering 3D phase-space trajectories...\n');
osc_colors = lines(num_oscillators);
fig_3d = figure('Name', '3D Phase Space');
ax3 = axes(fig_3d);
hold(ax3, 'on');
osc_labels_3d = cell(num_oscillators, 1);

for k = 1:num_oscillators
    plot3(ax3, x_states(:, k), y_states(:, k), z_states(:, k), ...
        'LineWidth', 1, 'Color', osc_colors(k, :));
    osc_labels_3d{k} = sprintf('Node %d', k);
end
grid(ax3, 'on');
view(ax3, 3);

set(ax3, 'FontSize', font_ticks);
xlabel(ax3, '$x(t)$', 'Interpreter', 'latex', 'FontSize', font_xlabel);
ylabel(ax3, '$y(t)$', 'Interpreter', 'latex', 'FontSize', font_ylabel);
zlabel(ax3, '$z(t)$', 'Interpreter', 'latex', 'FontSize', font_zlabel);
title(ax3, sprintf('3D Phase-Space Trajectories ($N = %d$, $k_s = %.2f$)', num_oscillators, ks_val), ...
    'Interpreter', 'latex', 'FontSize', font_title);
legend(ax3, osc_labels_3d, 'Interpreter', 'latex', 'FontSize', font_legend, 'Location', 'best');

% =========================================================================
% STEP 10: AUTOMATED EPS FIGURE EXPORT
% =========================================================================
fprintf('Saving figures to EPS format...\n');

folder_name = 'eps';
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

file_tag = sprintf('%s_%s_ks%.2f', FMtype, edges_input, ks_val);

filename_fig1 = fullfile(folder_name, sprintf('Motif_Trajectories_Error_%s.eps', file_tag));
exportgraphics(fig_combined, filename_fig1, 'ContentType', 'vector');
fprintf('  Saved: %s\n', filename_fig1);

filename_fig2 = fullfile(folder_name, sprintf('Motif_Manifold_%s.eps', file_tag));
exportgraphics(fig_planes, filename_fig2, 'ContentType', 'vector');
fprintf('  Saved: %s\n', filename_fig2);

filename_fig3 = fullfile(folder_name, sprintf('Motif_3DPhaseSpace_%s.eps', file_tag));
exportgraphics(fig_3d, filename_fig3, 'ContentType', 'vector');
fprintf('  Saved: %s\n', filename_fig3);

fprintf('Execution and EPS export finished successfully!\n');

% =========================================================================
% LOCAL EVENT FUNCTION FOR ODE45 BLOW-UP TERMINATION
% =========================================================================
function [value, isterminal, direction] = blowupEvent(~, Y, threshold)
    if any(isnan(Y)) || any(isinf(Y))
        value = 1;          % Positive: crossing from below is always detected
    else
        value = max(abs(Y)) - threshold;
    end
    isterminal = 1; 
    direction  = 1; 
end
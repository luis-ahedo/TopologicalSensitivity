% =========================================================================
% SCRIPT: MotifLLE_vs_ks_AllSchemes.m
% DESCRIPTION: For one topology, computes the largest transverse Lyapunov
%              exponent LLE(ks) for ks in I_ks = [0, kf] and every coupling
%              scheme CasoBs = 1 to 9, and plots the nine curves in a 3x3
%              tiledlayout (titles b_i -> c_j, shared axis labels).
%              For each scheme the variational equations on e = 0 are
%              rebuilt and odeMotifVariational.m is regenerated; the
%              exponents come from MaxTransLyaExp.
% =========================================================================
clear; clc; close all;

% =========================================================================
% STEP 1: CONFIGURATION
% =========================================================================
OsID = 1;                       % 1=Rossler, 2=Lorenz
num_oscillators = 3;
topology_idx = 1;               % Motif 1..13

kf     = 25;                    % I_ks = [0, kf]
N_ks   = 101;                   % Number of ks points
ks_vec = linspace(0, kf, N_ks);

opts = struct();
opts.tf          = 2000;        % Maximum averaging time
opts.dt          = 0.005;       % RK4 fixed step
opts.T_transient = 100;         % Transient of node 1
opts.num_LE      = 1;           % Largest exponent only
opts.seed        = 42;          % Same perturbation for every ks and scheme
opts.early_exit  = true;        % Stop once LLE has converged

use_parallel = false;           % true: parfor over ks (Parallel Computing Toolbox)

font_title     = 20;   % Global title (tiledlayout)
font_subtitle  = 20;   % Title of each subplot
font_xlabel    = 20;
font_ylabel    = 20;
font_ticks     = 15;

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
edges_input = topology_list{topology_idx};
edges       = eval(['edges_' edges_input]);

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
% STEP 3: SWEEP OVER COUPLING SCHEMES AND ks
% =========================================================================
LLE_all = NaN(9, N_ks);         % LLE_all(CasoBs, ks index)
scheme_title = strings(1, 9);
t_start = tic;

for CasoBs = 1:9
    [bs, cs] = SISObc(CasoBs, n);
    scheme_title(CasoBs) = sprintf('$\\mathbf{b}_%d \\to \\mathbf{c}_%d$', find(bs), find(cs));
    fprintf('\nCasoBs = %d / 9 (b_%d -> c_%d): building variational equations...\n', ...
        CasoBs, find(bs), find(cs));

    % Variational equations on e = 0 for this scheme -> odeMotifVariational.m
    buildVariationalODE(F, X, n, edges, bs, cs, Par);
    clear odeMotifVariational;
    rehash;
    if use_parallel
        pctRunOnAll clear odeMotifVariational;      % Workers pick up the new file
        pctRunOnAll rehash;
    end

    LLE_row = NaN(1, N_ks);
    if use_parallel
        parfor k_idx = 1:N_ks
            LLE_row(k_idx) = MaxTransLyaExp([OscPar, ks_vec(k_idx)], dim_Y, n, opts);
        end
    else
        for k_idx = 1:N_ks
            LLE_row(k_idx) = MaxTransLyaExp([OscPar, ks_vec(k_idx)], dim_Y, n, opts);
            fprintf('  ks = %6.3f  ->  LLE = %+.5f\n', ks_vec(k_idx), LLE_row(k_idx));
        end
    end
    LLE_all(CasoBs, :) = LLE_row;
end
fprintf('\nSweep finished in %.1f min\n', toc(t_start)/60);

% =========================================================================
% STEP 4: 3x3 TILED LAYOUT OF LLE(ks)
% =========================================================================
fig_lle = figure('Name', sprintf('LLE vs ks (%s - %s)', FMtype, edges_input));
tl = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for CasoBs = 1:9
    ax = nexttile(tl);
    plot(ax, ks_vec, LLE_all(CasoBs, :), 'b-', 'LineWidth', 1.5);
    hold(ax, 'on');
    yline(ax, 0, 'k--', 'LineWidth', 1.0);
    grid(ax, 'on');
    xlim(ax, [0 kf]);
    set(ax, 'FontSize', font_ticks);
    title(ax, scheme_title(CasoBs), 'Interpreter', 'latex', 'FontSize', font_subtitle);
end

% Shared axis labels for all subplots
xlabel(tl, '$k_s$', 'Interpreter', 'latex', 'FontSize', font_xlabel);
ylabel(tl, '$\lambda_{\max}$', 'Interpreter', 'latex', 'FontSize', font_ylabel);
title(tl, sprintf('%s - %s', osc_names{OsID}, motif_names{topology_idx}), ...
    'Interpreter', 'latex', 'FontSize', font_title);

% =========================================================================
% STEP 5: EPS EXPORT
% =========================================================================
folder_name = 'eps';
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end
filename_eps = fullfile(folder_name, sprintf('Motif_LLE_vs_ks_AllSchemes_%s_%s_kf%g.eps', ...
    FMtype, edges_input, kf));
exportgraphics(fig_lle, filename_eps, 'ContentType', 'vector');
fprintf('Saved: %s\n', filename_eps);

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

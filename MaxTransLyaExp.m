function [maxLyaExp, LyaExp, t_run, le_run] = MaxTransLyaExp(Parval, dim_Y, n, opts)
% =========================================================================
% FUNCTION: MaxTransLyaExp
% DESCRIPTION: Transverse Lyapunov exponents of the motif error system,
%              linearized on the synchronization manifold (e = 0).
%              Adapted from MaxLyaExpMotifsV2.m.
%
%   Integrates Z = [x1; delta e] with the generated odeMotifVariational.m
%   (built by MotifLLE_vs_ks_*.m): dx1/dt = F(x1), d(delta e)/dt =
%   J_E0(x1) * delta e. Fixed-step RK4; every column of the state matrix
%   is [x1; q_i], and QR re-orthonormalizes the q_i at every step.
%
% INPUTS
%   Parval : [a, b, c, ks]
%   dim_Y  : number of independent error states (6 = [de12; de13])
%   n      : states per node (3)
%   opts   : optional struct, fields (defaults):
%            tf (2000), dt (0.005), T_transient (100), num_LE (dim_Y),
%            x0 (rand(n,1)), seed ([] = no reseed), early_exit (true),
%            save_every (1, time units between stored running estimates),
%            odefun (@odeMotifVariational)
%
% OUTPUTS
%   maxLyaExp : largest transverse Lyapunov exponent (NaN on blow-up)
%   LyaExp    : exponents sorted in descending order (num_LE x 1)
%   t_run     : times of the stored running estimates
%   le_run    : running estimates (num_LE x numel(t_run))
% =========================================================================

% --- 1. Options ---
if nargin < 4, opts = struct(); end
tf          = getOpt(opts, 'tf',          2000);
dt          = getOpt(opts, 'dt',          0.005);
T_transient = getOpt(opts, 'T_transient', 100);
num_LE      = getOpt(opts, 'num_LE',      dim_Y);
seed        = getOpt(opts, 'seed',        []);
early_exit  = getOpt(opts, 'early_exit',  true);
save_every  = getOpt(opts, 'save_every',  1);
odefun      = getOpt(opts, 'odefun',      @odeMotifVariational);
if ~isempty(seed), rng(seed, 'twister'); end   % Explicit generator: parallel workers default
                                               % to another one, so rng(seed) alone would give
                                               % other random numbers on a worker than on the client
x0          = getOpt(opts, 'x0',          rand(n, 1));

hdt       = 0.5 * dt;
ns        = round(tf / dt);                     % Number of RK4 steps
save_k    = max(1, round(save_every / dt));     % Steps between stored estimates
maxLyaExp = NaN;
LyaExp    = NaN(num_LE, 1);
t_run     = zeros(1, floor(ns / save_k));
le_run    = NaN(num_LE, floor(ns / save_k));
i_save    = 0;

f = @(tk, Z) odefun(tk, Z, Parval);

% --- 2. Transient of node 1 (perturbations set to zero) ---
Z = [x0; zeros(dim_Y, 1)];
for k = 1:round(T_transient / dt)
    Z = rk4Step(f, 0, Z, dt, hdt);
end
if any(~isfinite(Z))
    return;
end
x1 = Z(1:n);

% --- 3. Initial state matrix [x1 ... x1; Q] ---
% First perturbation consistent with the loop closure de23 = de13 - de12
Y0 = randn(dim_Y, num_LE);
if dim_Y == 3*n
    Y0(2*n+1:3*n, 1) = Y0(n+1:2*n, 1) - Y0(1:n, 1);
end
[Q, ~] = qr(Y0, 0);
xk     = [repmat(x1, 1, num_LE); Q];
LyaEx  = zeros(num_LE, 1);

% --- 4. Convergence test configuration ---
check_interval = round(25 / dt);                % Check every 25 time units
min_time_steps = round(0.2 * ns);               % No exit before 20% of tf
window_history = 10;
conv_tol       = 1e-4;
max_le_history = [];
T_used         = tf;

% --- 5. Integration loop ---
for k = 1:ns
    tk = (k - 1) * dt;
    xk = rk4Step(f, tk, xk, dt, hdt);

    % Down-sampled blow-up safeguard
    if mod(k, 500) == 0
        if any(~isfinite(xk(:))) || any(abs(xk(:)) > 1e10)
            return;
        end
    end

    % QR (Gram-Schmidt) re-orthonormalization
    [Q, R] = qr(xk(n+1:end, :), 0);
    LyaEx  = LyaEx + log(abs(diag(R)));
    xk(n+1:end, :) = Q;

    % Running estimates (for convergence plots)
    if mod(k, save_k) == 0
        i_save = i_save + 1;
        t_run(i_save)     = k * dt;
        le_run(:, i_save) = LyaEx / (k * dt);
    end

    % Dual-condition early termination
    if early_exit && k > min_time_steps && mod(k, check_interval) == 0
        current_time   = k * dt;
        current_max_le = max(LyaEx / current_time);
        max_le_history = [max_le_history, current_max_le]; %#ok<AGROW>
        if numel(max_le_history) > window_history
            max_le_history(1) = [];
        end
        if numel(max_le_history) == window_history
            % Converged (stable or unstable) or clearly chaotic transverse dynamics
            if std(max_le_history) < conv_tol * (abs(current_max_le) + 1) || ...
               all(max_le_history > 0.01)
                T_used = current_time;
                break;
            end
        end
    end
end

% --- 6. Final normalization ---
LyaExp    = sort(LyaEx / T_used, 'descend');
maxLyaExp = LyaExp(1);
t_run     = t_run(1:i_save);
le_run    = le_run(:, 1:i_save);
end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function Z = rk4Step(f, tk, Z, dt, hdt)
    k1 = f(tk,       Z);
    k2 = f(tk + hdt, Z + hdt * k1);
    k3 = f(tk + hdt, Z + hdt * k2);
    k4 = f(tk + dt,  Z + dt  * k3);
    Z  = Z + (k1 + 2*k2 + 2*k3 + k4) * (dt / 6);
end

function v = getOpt(opts, name, default)
    if isfield(opts, name) && ~isempty(opts.(name))
        v = opts.(name);
    else
        v = default;
    end
end

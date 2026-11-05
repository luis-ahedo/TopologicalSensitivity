function [x_ref, x_ref_std, X_rms] = AttractorScale(OsID, recompute)
% =========================================================================
% FUNCTION: AttractorScale.m
% DESCRIPTION: Attractor scale x_ref of Eq. (eq_normalized_sensitivity),
%              one value per oscillator:
%                  x_ref = RMS(x) = ( sum_l rms(x_l)^2 )^(1/2)
%              (vector RMS, Eq. eq:rms-vector) of the nominal node
%              trajectory over the steady-state window [t_transient, T_end]
%              of the RMS analysis, averaged over n_ic initial conditions
%              uniform in [0, 1]^3. The normalized sensitivity is
%                  sbar_ij = s_ij * q_j / x_ref.
%              The value is stored in mat/AttractorScale_<Osc>.mat and read
%              from there in later calls (recompute = true runs it again).
%              The random-number state of the caller is restored.
% OUTPUT:      x_ref     : mean of the n_ic values
%              x_ref_std : standard deviation of the n_ic values
%              X_rms     : the n_ic values
% =========================================================================
if nargin < 2, recompute = false; end

osc_list = {'Rossler', 'Lorenz'};
FMtype   = osc_list{OsID};
switch FMtype
    case 'Rossler'
        OscPar = [0.2, 0.2, 5.7];
        Fnode  = @(t, x) [-x(2) - x(3); x(1) + OscPar(1)*x(2); OscPar(2) + x(3)*(x(1) - OscPar(3))];
    case 'Lorenz'
        OscPar = [10, 28, 2];     % same values as OscillatorType.m
        Fnode  = @(t, x) [OscPar(1)*(x(2) - x(1)); x(1)*(OscPar(2) - x(3)) - x(2); x(1)*x(2) - OscPar(3)*x(3)];
end

% Same window, step and solver settings as the RMS analysis (Appendix A.5)
t_transient = 500;
T_end       = 2500;
dt          = 0.01;
n_ic        = 10;
base_seed   = 42;
ode_RelTol  = 1e-8;
ode_AbsTol  = 1e-10;

cache_file = fullfile('mat', sprintf('AttractorScale_%s.mat', FMtype));
if ~recompute && isfile(cache_file)
    L = load(cache_file);
    if isequal(L.OscPar, OscPar) && L.t_transient == t_transient && L.T_end == T_end ...
            && L.dt == dt && L.n_ic == n_ic
        x_ref = L.x_ref;  x_ref_std = L.x_ref_std;  X_rms = L.X_rms;
        return;
    end
end

rng_state = rng;                                    % restored at the end
rng(base_seed, 'twister');
t_samp = t_transient + (1:round((T_end - t_transient) / dt)) * dt;
opts   = odeset('RelTol', ode_RelTol, 'AbsTol', ode_AbsTol);
X_rms  = NaN(n_ic, 1);
for r = 1:n_ic
    x0 = rand(3, 1);
    [~, X] = ode45(Fnode, [0, t_samp], x0, opts);
    Xs = X(2:end, :);                               % samples in the window
    X_rms(r) = sqrt(sum(mean(Xs.^2, 1)));           % vector RMS
end
rng(rng_state);
x_ref     = mean(X_rms);
x_ref_std = std(X_rms);

if ~exist('mat', 'dir'), mkdir('mat'); end
save(cache_file, 'x_ref', 'x_ref_std', 'X_rms', 'OscPar', 'FMtype', 't_transient', 'T_end', ...
    'dt', 'n_ic', 'base_seed', 'ode_RelTol', 'ode_AbsTol');
fprintf('AttractorScale (%s): x_ref = %.4f +/- %.4f (%d initial conditions), -20 log10(x_ref) = %.2f dB\n', ...
    FMtype, x_ref, x_ref_std, n_ic, -20*log10(x_ref));
end

% =========================================================================
% SCRIPT: SensitivityNorm_Consolidate.m
% DESCRIPTION: Joins the per-combination files written by
%              SensitivityNorm_vs_ks_AllTopologies.m (or SensitivityNorm_vs_ks.m),
%                  mat/SensitivityNorm_vs_ks_<Osc>_<topology>_B<C>_N<N>_T<T_end>.mat
%              (13 topologies x 9 schemes), into one unified file,
%                  mat/SensitivityNorm_vs_ks_<Osc>_AllTopologies_N<N>_T<T_end>.mat
%              with everything needed to redraw the figures of
%              SensitivityNorm_vs_ks.m for any topology and scheme:
%                SN(topology, CasoBs)  13 x 9 struct, same field names as
%                  the single-combination script (ks_sweep, lle_sweep,
%                  norm_dB, norm_unc, ks_ok, no_influence, S_worst, S_mean,
%                  ks_best, S_best, k_best, ks_plot, ks_a, ks_b, lam_min, ...)
%                  plus lle_row (lambda_max over the LLE grid, panel (a))
%                  and status: 'ok', 'all_flagged', 'no_interval' (no valid
%                  interval in the LLE file, or removed by the RMS file: no file
%                  exists), 'missing' (valid case whose file is not there yet);
%                ks_lle, kf, LLE_array   LLE grid for panel (a);
%                summary tables 13 x 9  (status_tab, ks_a_tab, ks_b_tab,
%                  ks_best_tab, S_best_tab, n_notconv_tab, n_flagged_tab,
%                  runs_tab, time_tab);
%                names and labels (topology_list, motif_names, scheme_bc,
%                  osc_names, par_tex, err_tex) and the run settings.
%              Topological sensitivity index of each mismatch q (status 'ok'):
%                Step 1: S_q(e)(ks) = sqrt( (1/3) sum_ij S_q(e_ij)^2 )
%                        (linear values, no influence = 0, fixed 3)
%                Step 2: <S_q> = mean over the valid ks of 20 log10 S_q(e)(ks)
%                        (valid ks: RMS check passed AND convergence test met,
%                        the same rule as the TopologicalSensitivity_* figures)
%              Per combination: SN.TS_dB (1 x 6), SN.TS_curve_dB (N x 6,
%              step-1 curve in dB), SN.TS_unc (1 x 6, 1.2533 * std / sqrt(n) of
%              the index per run, runs valid at every valid ks), SN.TS_gap
%              (|TS_dB - median of the per-run indices|), SN.TS_Nr (runs);
%              tables 13 x 9 x 6: TS_tab, TS_unc_tab, TS_gap_tab.
%              The per-run RMS values (RMS_all, rms_e) stay in the
%              per-combination files unless include_raw = true.
%              Run it on its own (set OsID, N_sweep, time_profile below), or
%              it is called at the end of SensitivityNorm_vs_ks_AllTopologies.m
%              (consolidate_at_end), which passes its settings.
% =========================================================================
if exist('consolidate_OsID', 'var')
    OsID    = consolidate_OsID;             % Called from SensitivityNorm_vs_ks_AllTopologies.m
    N_sweep = consolidate_N_sweep;
    T_end   = consolidate_T_end;
else
    clear; clc; close all;
    % =====================================================================
    % CONFIGURATION (stand-alone use)
    % =====================================================================
    OsID         = 2;                       % 1=Rossler, 2=Lorenz
    N_sweep      = 25;                      % Same N_sweep as the computation
    time_profile = 'medium';                % Same time_profile as the computation
    switch lower(time_profile)
        case 'short',  T_end = 1000;
        case 'medium', T_end = 2500;
        case 'long',   T_end = 5000;
        otherwise, error('time_profile must be ''short'', ''medium'' or ''long''.');
    end
end
% Rule for a valid case: the same as SensitivityNorm_vs_ks_AllTopologies.m
% (passed by it in consolidate_rule; defaults below for stand-alone use).
% Only valid cases have a file; a valid case without a file is 'missing'.
if exist('consolidate_rule', 'var')
    rule = consolidate_rule;
else
    rule = struct('interval_pick', 0, 'n_decay', 5, 'min_width', 0.1, ...
                  'use_rms_file', true, 'rms_file_tol', []);
end
t_transient_rule = T_end / 5;               % every time profile uses t_transient = T_end / 5
lam_min_rule     = rule.n_decay / t_transient_rule;
allow_incomplete = false;                   % true: save even if some combinations are missing
include_raw      = false;                   % true: also copy RMS_all and rms_e (large)
mat_folder       = 'mat';

OscType = {'Rossler', 'Lorenz', 'Chen', 'HR', 'Duffing', 'vanderPol'};
FMtype  = OscType{OsID};                    % Same names as OscillatorType.m

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
num_topo   = numel(topology_list);
num_scheme = 9;
osc_names  = {'R\"{o}ssler', 'Lorenz'};
par_tex    = {'\delta a_2', '\delta b_2', '\delta c_2', '\delta a_3', '\delta b_3', '\delta c_3'};
err_tex    = {'12', '13', '23'};

% =========================================================================
% LLE GRID (panel (a) of the robustness figure)
% =========================================================================
lle_file = fullfile(mat_folder, sprintf('LLE_vs_ks_%s_AllTopologies.mat', FMtype));
if ~isfile(lle_file)
    error('%s not found.', lle_file);
end
L = load(lle_file, 'LLE_array', 'ks_vec', 'kf');
LLE_array = L.LLE_array;
ks_lle    = L.ks_vec(:).';
kf        = L.kf;

% =========================================================================
% STRUCT TEMPLATE (same field names as SensitivityNorm_vs_ks.m)
% =========================================================================
data_fields = {'ks_sweep', 'lle_sweep', 'norm_dB', 'norm_unc', 'norm_dB_all', 'norm_unc_all', ...
    'comp_dB', 'ks_ok', 'no_influence', 'S_worst', 'S_mean', 'ks_best', 'S_best', 'k_best', ...
    'ks_plot', 'ks_a', 'ks_b', 'I_all', 'r_pick', 'lam_min', 'n_used', 'conv_ks', 'n_blown', ...
    'rms_e_med', 'run_time_min', 'date_computed', 'reason'};
ts_fields   = {'TS_dB', 'TS_curve_dB', 'TS_unc', 'TS_gap', 'TS_Nr'};   % topological sensitivity index
raw_fields  = {'RMS_all', 'rms_e'};
tmpl = struct('status', 'missing', 'topology_idx', NaN, 'CasoBs', NaN, 'scheme_bc', [NaN NaN], ...
    'motif_name', '', 'file', '', 'lle_row', []);
for f = data_fields, tmpl.(f{1}) = []; end
for f = ts_fields,   tmpl.(f{1}) = []; end
if include_raw
    for f = raw_fields, tmpl.(f{1}) = []; end
end
SN = repmat(tmpl, num_topo, num_scheme);
xref = NaN;                                   % attractor scale x_ref, read from the files

% Settings checked for consistency between files
set_keys = {'max_runs', 'batch_size', 'min_runs', 'conv_window', 'conv_tol_dB', 'conv_sem_dB', ...
    't_transient', 'T_end', 'dt_sample', 'base_seed', 'n_decay', 'use_rms_file', 'rms_file_tol', ...
    'use_internal_rms_check', 'rms_e_tol', 'time_profile'};
settings = struct();
set_from = '';

% =========================================================================
% READ THE PER-COMBINATION FILES
% =========================================================================
status_tab    = strings(num_topo, num_scheme);
ks_a_tab      = NaN(num_topo, num_scheme);
ks_b_tab      = NaN(num_topo, num_scheme);
ks_best_tab   = NaN(num_topo, num_scheme);
S_best_tab    = NaN(num_topo, num_scheme);
n_notconv_tab = NaN(num_topo, num_scheme);
n_flagged_tab = NaN(num_topo, num_scheme);
runs_tab      = NaN(num_topo, num_scheme);
time_tab      = NaN(num_topo, num_scheme);
TS_tab        = NaN(num_topo, num_scheme, 6);   % <S_q> (dB), q = da2 db2 dc2 da3 db3 dc3
TS_unc_tab    = NaN(num_topo, num_scheme, 6);
TS_gap_tab    = NaN(num_topo, num_scheme, 6);
for ti = 1:num_topo
    Rr = [];
    if rule.use_rms_file
        rms_file_r = fullfile(mat_folder, sprintf('RMS_vs_ks_%s_%s_T%d.mat', FMtype, topology_list{ti}, T_end));
        if isfile(rms_file_r)
            Rr = load(rms_file_r, 'RMS_topo', 'ks_vec', 'sync_tol');
            rms_tol_r = rule.rms_file_tol;
            if isempty(rms_tol_r), rms_tol_r = Rr.sync_tol; end
        else
            warning('%s not found: validity of %s decided from the LLE file only.', rms_file_r, topology_list{ti});
        end
    end
    for CasoBs = 1:num_scheme
        s = tmpl;
        s.topology_idx = ti;
        s.CasoBs       = CasoBs;
        s.motif_name   = motif_names{ti};
        s.lle_row      = reshape(LLE_array(ti, CasoBs, :), 1, []);
        fname = fullfile(mat_folder, sprintf('SensitivityNorm_vs_ks_%s_%s_B%d_N%d_T%d.mat', ...
            FMtype, topology_list{ti}, CasoBs, N_sweep, T_end));
        if ~isfile(fname)
            % No file: a case without a valid interval (not computed) or a missing one
            LI = lleInterval(ks_lle, s.lle_row, rule.interval_pick, lam_min_rule, rule.min_width, false);
            reason_r = LI.reason;
            if isempty(reason_r) && ~isempty(Rr)
                [~, ~, reason_r] = rmsTrim(LI.ks_a, LI.ks_b, Rr.ks_vec(:).', Rr.RMS_topo(CasoBs, :), ...
                    rms_tol_r, rule.min_width, false);
            end
            if ~isempty(reason_r)
                s.status = 'no_interval';
                s.reason = reason_r;
            end
        end
        if isfile(fname)
            T = load(fname);
            s.file   = fname;
            s.status = 'ok';                        % files of SensitivityNorm_vs_ks.m have no status
            if isfield(T, 'status'), s.status = T.status; end
            if strcmp(s.status, 'ok')               % same attractor scale x_ref in every file
                if ~isfield(T, 'xref')
                    error('%s has no xref field (run ShiftSensitivity_xref.m first).', fname);
                elseif isnan(xref)
                    xref = T.xref;
                elseif abs(T.xref - xref) > 1e-12 * xref
                    error('%s: xref = %.6g differs from %.6g of the other files.', fname, T.xref, xref);
                end
            end
            if isfield(T, 'scheme_bc'), s.scheme_bc = T.scheme_bc; end
            for f = data_fields
                if isfield(T, f{1}), s.(f{1}) = T.(f{1}); end
            end
            if include_raw
                for f = raw_fields
                    if isfield(T, f{1}), s.(f{1}) = T.(f{1}); end
                end
            end
            % Fields not saved by SensitivityNorm_vs_ks.m
            if strcmp(s.status, 'ok')
                if isempty(s.k_best)
                    s.k_best = find(s.ks_sweep == s.ks_best, 1);
                end
                if isempty(s.ks_plot)
                    s.ks_plot = [min(s.ks_sweep(s.ks_ok)), max(s.ks_sweep(s.ks_ok))];
                end
                % Topological sensitivity index (needs RMS_all for the uncertainty)
                RMS_raw = [];
                if isfield(T, 'RMS_all'), RMS_raw = T.RMS_all; end
                % Valid ks: RMS check passed AND convergence test met (same rule as the figures)
                ks_valid = s.ks_ok(:);
                if ~isempty(s.conv_ks), ks_valid = ks_valid & s.conv_ks(:); end
                TS = topoSensIndex(s.norm_dB, s.no_influence, ks_valid, RMS_raw);
                for f = ts_fields, s.(f{1}) = TS.(f{1}); end
                TS_tab(ti, CasoBs, :)     = TS.TS_dB;
                TS_unc_tab(ti, CasoBs, :) = TS.TS_unc;
                TS_gap_tab(ti, CasoBs, :) = TS.TS_gap;
            end
            % Settings consistency
            if isempty(set_from)
                for f = set_keys
                    if isfield(T, f{1}), settings.(f{1}) = T.(f{1}); end
                end
                set_from = fname;
            else
                diff_keys = {};
                for f = set_keys
                    if isfield(T, f{1}) && isfield(settings, f{1}) && ~isequal(T.(f{1}), settings.(f{1}))
                        diff_keys{end+1} = f{1};    %#ok<AGROW>
                    end
                end
                if ~isempty(diff_keys)
                    warning('%s: %s differ from %s.', fname, strjoin(diff_keys, ', '), set_from);
                end
            end
            % Summary tables
            ks_a_tab(ti, CasoBs) = scalarOrNaN(s.ks_a);
            ks_b_tab(ti, CasoBs) = scalarOrNaN(s.ks_b);
            if ~isempty(s.run_time_min), time_tab(ti, CasoBs) = s.run_time_min; end
            if any(strcmp(s.status, {'ok', 'all_flagged'}))
                ks_best_tab(ti, CasoBs)   = s.ks_best;
                S_best_tab(ti, CasoBs)    = s.S_best;
                n_notconv_tab(ti, CasoBs) = nnz(~s.conv_ks);
                n_flagged_tab(ti, CasoBs) = nnz(~s.ks_ok);
                runs_tab(ti, CasoBs)      = sum(s.n_used);
            end
        end
        if any(isnan(s.scheme_bc))
            [bs, cs] = SISObc(CasoBs, 3);
            s.scheme_bc = [find(bs), find(cs)];
        end
        SN(ti, CasoBs)        = s;
        status_tab(ti, CasoBs) = s.status;
    end
end
scheme_bc = reshape([SN(1, :).scheme_bc], 2, []).';
scheme_names = compose('b%d_to_c%d', scheme_bc);

% =========================================================================
% SUMMARY
% =========================================================================
n_missing = nnz(status_tab == "missing");
fprintf('\n%s, N_sweep = %d, T_end = %d: %d ok, %d all flagged, %d no interval, %d missing (of %d)\n', ...
    FMtype, N_sweep, T_end, nnz(status_tab == "ok"), nnz(status_tab == "all_flagged"), ...
    nnz(status_tab == "no_interval"), n_missing, num_topo * num_scheme);
fprintf('\n%-32s %-9s %-18s %9s %10s %8s %8s %7s\n', 'Topology', 'Scheme', 'Swept ks', 'ks_best', ...
    'S_best dB', 'not conv', 'flagged', 'min');
for ti = 1:num_topo
    for CasoBs = 1:num_scheme
        s = SN(ti, CasoBs);
        sch = sprintf('b%d->c%d', s.scheme_bc(1), s.scheme_bc(2));
        switch s.status
            case {'ok', 'all_flagged'}
                fprintf('%-32s %-9s [%6.3f, %6.3f] %9.3f %10.2f %8d %8d %7.1f\n', motif_names{ti}, sch, ...
                    s.ks_a, s.ks_b, s.ks_best, s.S_best, n_notconv_tab(ti, CasoBs), ...
                    n_flagged_tab(ti, CasoBs), time_tab(ti, CasoBs));
            case 'no_interval'
                fprintf('%-32s %-9s no sweep: %s\n', motif_names{ti}, sch, s.reason);
            otherwise
                fprintf('%-32s %-9s missing\n', motif_names{ti}, sch);
        end
    end
end

% Topological sensitivity index <S_q> (dB) +/- uncertainty
fprintf('\nTopological sensitivity index <S_q> (dB, mean over valid ks of the fixed-3 combination) +/- uncertainty\n');
fprintf('%-32s %-9s', 'Topology', 'Scheme');
fprintf(' %14s', 'da2', 'db2', 'dc2', 'da3', 'db3', 'dc3');
fprintf(' %8s\n', 'max gap');
for ti = 1:num_topo
    for CasoBs = 1:num_scheme
        s = SN(ti, CasoBs);
        if isempty(s.TS_dB), continue; end
        fprintf('%-32s %-9s', motif_names{ti}, sprintf('b%d->c%d', s.scheme_bc(1), s.scheme_bc(2)));
        fprintf(' %7.2f +/-%4.2f', [s.TS_dB; s.TS_unc]);
        fprintf(' %8.2f\n', max(s.TS_gap));
    end
end

% =========================================================================
% CHECK AND SAVE
% =========================================================================
if n_missing > 0
    [ti_m, b_m] = find(status_tab == "missing");
    missing_topo = unique(ti_m).';
    fprintf('\nMissing combinations in topologies: %s\n', strjoin(topology_list(missing_topo), ', '));
    if ~allow_incomplete
        warning('%d of %d combinations missing: unified file NOT written (set allow_incomplete = true to write it anyway).', ...
            n_missing, num_topo * num_scheme);
        return;
    end
end

filename_all = fullfile(mat_folder, sprintf('SensitivityNorm_vs_ks_%s_AllTopologies_N%d_T%d.mat', FMtype, N_sweep, T_end));
date_consolidated = datestr(now, 'yyyy-mm-dd HH:MM');
mat_version = '-v7';
if include_raw, mat_version = '-v7.3'; end   % raw runs can exceed 2 GB
save(filename_all, 'SN', 'status_tab', 'ks_a_tab', 'ks_b_tab', 'ks_best_tab', 'S_best_tab', ...
    'n_notconv_tab', 'n_flagged_tab', 'runs_tab', 'time_tab', 'TS_tab', 'TS_unc_tab', 'TS_gap_tab', 'LLE_array', 'ks_lle', 'kf', ...
    'topology_list', 'motif_names', 'scheme_bc', 'scheme_names', 'osc_names', 'par_tex', 'err_tex', ...
    'settings', 'FMtype', 'OsID', 'N_sweep', 'T_end', 'lle_file', 'include_raw', 'date_consolidated', 'xref', mat_version);
fprintf('\nSaved: %s (total computing time %.1f h)\n', filename_all, sum(time_tab(:), 'omitnan') / 60);

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
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

function [i1, i2] = longestRun(mask)
    % First and last index of the longest run of true values in mask
    d  = diff([0, mask(:).', 0]);
    on = find(d == 1);
    of = find(d == -1) - 1;
    [~, r] = max(of - on);
    i1 = on(r);
    i2 = of(r);
end

function v = scalarOrNaN(x)
    if isempty(x), v = NaN; else, v = x; end
end

% -------------------------------------------------------------------------
% Topological sensitivity index (same computation as TopologicalSensitivity_Test.m)
%   norm_dB (N x 18, v = (j-1)*3 + p), no_influence (1 x 18), ks_ok (N x 1: valid ks),
%   RMS_all (N x max_runs x 54) or [] (then no uncertainty).
% -------------------------------------------------------------------------
function TS = topoSensIndex(norm_dB, no_influence, ks_ok, RMS_all)
    Nq = 6;  Np = 3;  n = 3;  n_err = n * Np;
    N    = size(norm_dB, 1);
    ok   = logical(ks_ok(:));
    % A vector that is NaN at a valid ks while other vectors are finite has
    % all its coefficients structurally zero there (STEP D of the profile):
    % it counts as 0, as no_influence does at every ks. A ks where every
    % vector is NaN (no valid runs) is left out.
    ok   = ok & any(isfinite(norm_dB), 2);
    kidx = find(ok);
    NI   = reshape(no_influence, Np, Nq).';            % NI(j, p)
    TS   = struct('TS_dB', NaN(1, Nq), 'TS_curve_dB', NaN(N, Nq), 'TS_unc', NaN(1, Nq), ...
                  'TS_gap', NaN(1, Nq), 'TS_Nr', 0);
    if isempty(kidx), return; end

    % Step 1 (linear, no influence = 0, fixed 3) and step 2 (mean in dB)
    S = NaN(N, Nq, Np);
    for j = 1:Nq
        for p = 1:Np
            if NI(j, p)
                S(:, j, p) = 0;
            else
                S(:, j, p) = 10.^(norm_dB(:, (j-1)*Np + p) / 20);
                S(ok & isnan(norm_dB(:, (j-1)*Np + p)), j, p) = 0;   % zero vector at this ks
            end
        end
    end
    S(~ok, :, :)   = NaN;
    S3             = sqrt(mean(S.^2, 3));              % N x Nq
    TS.TS_curve_dB = 20*log10(S3);
    TS.TS_dB       = mean(TS.TS_curve_dB(kidx, :), 1);

    % Uncertainty from the index per run (runs valid at every valid ks)
    if isempty(RMS_all), return; end
    max_runs = size(RMS_all, 2);
    Sr = NaN(N, max_runs, Nq, Np);
    for j = 1:Nq
        for p = 1:Np
            if NI(j, p)
                Sr(:, :, j, p) = 0;
            else
                cols = (j-1)*n_err + (p-1)*n + (1:n);
                Sr(:, :, j, p) = sqrt(sum(RMS_all(:, :, cols).^2, 3));
                Sr(ok & isnan(norm_dB(:, (j-1)*Np + p)), :, j, p) = 0;
            end
        end
    end
    S3r    = sqrt(mean(Sr.^2, 4));                     % N x max_runs x Nq
    Ir     = reshape(mean(20*log10(S3r(kidx, :, :)), 1), max_runs, Nq);
    run_ok = all(isfinite(Ir), 2);
    TS.TS_Nr = nnz(run_ok);
    if TS.TS_Nr < 2, return; end
    TS.TS_unc = 1.2533 * std(Ir(run_ok, :), 0, 1) / sqrt(TS.TS_Nr);
    TS.TS_gap = abs(TS.TS_dB - median(Ir(run_ok, :), 1));
end

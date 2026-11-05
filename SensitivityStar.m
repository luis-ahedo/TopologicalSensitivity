% =========================================================================
% SCRIPT: SensitivityStar.m
% DESCRIPTION: Global sensitivity index at maximal transverse stability,
%              S_q*, for one oscillator, from the unified file
%                  mat/SensitivityNorm_vs_ks_<Osc>_AllTopologies_N<N>_T<T_end>.mat
%              made by SensitivityNorm_Consolidate.m (no new simulations).
%              For every topology and coupling scheme with a result:
%                  ks*   = argmin over the sweep of lambda_max(ks)
%                  band  = { ks : lambda_max(ks) <= band_frac * lambda_max(ks*) }
%                  S_q*  = mean over the band of 20 log10 Sbar_q(e)(ks)
%              where Sbar_q(e) is the step-1 curve sqrt((1/3) sum_ij S_q(e_ij)^2)
%              of SensitivityNorm_Consolidate.m (SN.TS_curve_dB). Unlike
%              <S_q>, S_q* does not depend on the width of the
%              synchronization interval, on the excluded edges or on the
%              cutoff, as long as ks* lies inside the swept interval.
%              Cases whose ks* is the last swept point of an interval cut at
%              the end of the computed range (Rossler: kf = 5, Lorenz:
%              ks_max = 20) are flagged: there lambda_max still decreases at
%              the cutoff, so S_q* depends on it.
%              Output: command-window table,
%                  mat/SensitivityStar_<Osc>_N<N>_T<T_end>.mat
%                  tex/SensitivityStar_<Osc>_AllSchemes_N<N>_T<T_end>.tex
% =========================================================================
clear; clc; close all;

% =========================================================================
% STEP 1: CONFIGURATION
% =========================================================================
OsID         = 1;               % 1=Rossler, 2=Lorenz
N_sweep      = 25;              % Same N_sweep as the computation
time_profile = 'medium';        % Same time_profile as the computation: 'short' | 'medium' | 'long'
scheme_run   = [];              % [] = every scheme with a result; or CasoBs values, e.g. [1 5]
band_frac    = 0.9;             % band: lambda_max <= band_frac * min(lambda_max); 1 = the single point ks*
ks_cut_table = [5 20];          % per OsID: end of the computed ks range (Rossler kf, Lorenz ks_max)

mat_folder = 'mat';
tex_folder = 'tex';
write_table = true;
num_fmt     = '%.1f';           % Number format of S_q* (dB)
ks_fmt      = '%.2f';           % Number format of ks*

OscType   = {'Rossler', 'Lorenz', 'Chen', 'HR', 'Duffing', 'vanderPol'};   % Same names as OscillatorType.m
FMtype    = OscType{OsID};
osc_tex   = {'R\"{o}ssler', 'Lorenz', 'Chen', 'Hindmarsh-Rose', 'Duffing', 'van der Pol'};
topo_short = {'Lin. chain', 'Lin. aux.', 'Lin. relay', 'Lin. mix. mast.', 'Lin. comp.', ...
    'Lin. mix. slav.', 'Ring bidir. chain', 'Ring bidir. dbl. mix.', 'Ring bidir. mix.', ...
    'Ring unidir. chain', 'Ring bidir. aux.', 'Ring unidir. aux.', 'Ring bidir. comp.'};
par_txt     = {'da2', 'db2', 'dc2', 'da3', 'db3', 'dc3'};
par_tex_tab = {'\delta \mathcal{a}_2', '\delta \mathcal{b}_2', '\delta \mathcal{c}_2', ...
    '\delta \mathcal{a}_3', '\delta \mathcal{b}_3', '\delta \mathcal{c}_3'};
switch lower(time_profile)
    case 'short',  T_end = 1000;
    case 'medium', T_end = 2500;
    case 'long',   T_end = 5000;
    otherwise, error('time_profile must be ''short'', ''medium'' or ''long''.');
end
first_ring = 7;                 % Topologies 7..13 are rings
ks_cut     = ks_cut_table(OsID);

% =========================================================================
% STEP 2: LOAD THE UNIFIED FILE (made by SensitivityNorm_Consolidate.m)
% =========================================================================
filename_mat = fullfile(mat_folder, sprintf('SensitivityNorm_vs_ks_%s_AllTopologies_N%d_T%d.mat', FMtype, N_sweep, T_end));
if ~isfile(filename_mat)
    error('%s not found: run SensitivityNorm_Consolidate.m first.', filename_mat);
end
U = load(filename_mat, 'SN', 'TS_tab', 'scheme_bc', 'topology_list');
num_topo = numel(U.topology_list);
fprintf('Loaded: %s\n', filename_mat);

has_result = any(isfinite(U.TS_tab), 3);                   % 13 x 9, same cases as <S_q>
if isempty(scheme_run)
    scheme_run = find(any(has_result, 1));
end

% =========================================================================
% STEP 3: S_q* PER TOPOLOGY AND SCHEME
% =========================================================================
S_star   = NaN(num_topo, 9, 6);     % S_q* (dB)
ks_star  = NaN(num_topo, 9);        % argmin of lambda_max
lam_star = NaN(num_topo, 9);        % min of lambda_max
band_lo  = NaN(num_topo, 9);        % band limits and number of ks in it
band_hi  = NaN(num_topo, 9);
band_n   = zeros(num_topo, 9);
at_cut   = false(num_topo, 9);      % ks* is the last point of an interval cut at ks_cut
conv_cut = false(num_topo, 9);      % min of lambda_max at a ks without a converged index

for CasoBs = scheme_run
    for ti = 1:num_topo
        if ~has_result(ti, CasoBs), continue; end
        R   = U.SN(ti, CasoBs);
        ks  = R.ks_sweep(:);
        lam = R.lle_sweep(:);
        C   = R.TS_curve_dB;                                % N_sweep x 6 (dB)
        ok  = isfinite(lam) & any(isfinite(C), 2);
        if ~any(ok), continue; end
        lam_v = lam;
        lam_v(~ok) = NaN;                                   % min() ignores NaN
        [lam_min, i0] = min(lam_v);
        band = ok & (lam <= band_frac * lam_min);
        S_star(ti, CasoBs, :) = mean(C(band, :), 1, 'omitnan');
        ks_star(ti, CasoBs)   = ks(i0);
        lam_star(ti, CasoBs)  = lam_min;
        band_lo(ti, CasoBs)   = min(ks(band));
        band_hi(ti, CasoBs)   = max(ks(band));
        band_n(ti, CasoBs)    = nnz(band);
        at_cut(ti, CasoBs)    = (i0 == find(ok, 1, 'last')) && abs(ks(end) - ks_cut) < 1e-6;
        [~, i_all] = min(lam);                              % over every ks of the sweep
        conv_cut(ti, CasoBs)  = ~ok(i_all);
    end
end

% =========================================================================
% STEP 4: COMMAND WINDOW TABLE
% =========================================================================
for CasoBs = scheme_run
    bc = U.scheme_bc(CasoBs, :);
    fprintf('\nS_q* (dB), %s, b%d -> c%d (CasoBs = %d), band_frac = %.2f\n', FMtype, bc(1), bc(2), CasoBs, band_frac);
    fprintf('%-22s %6s %7s %-15s %3s', 'Topology', 'ks*', 'lam*', 'band', 'n');
    fprintf(' %7s', par_txt{:});
    fprintf('\n');
    for ti = 1:num_topo
        if ~has_result(ti, CasoBs), continue; end
        fprintf('%-22s %6.2f %7.3f %-15s %3d', topo_short{ti}, ks_star(ti, CasoBs), lam_star(ti, CasoBs), ...
            sprintf('[%.2f, %.2f]', band_lo(ti, CasoBs), band_hi(ti, CasoBs)), band_n(ti, CasoBs));
        fprintf(' %7.1f', S_star(ti, CasoBs, :));
        if at_cut(ti, CasoBs), fprintf('  (ks* at the cutoff %g)', ks_cut); end
        if conv_cut(ti, CasoBs), fprintf('  (min of lambda_max at a ks without a converged index)'); end
        fprintf('\n');
    end
end
fprintf('\n%d of %d cases have ks* at the cutoff.\n', nnz(at_cut(:, scheme_run)), nnz(has_result(:, scheme_run)));

filename_out = fullfile(mat_folder, sprintf('SensitivityStar_%s_N%d_T%d.mat', FMtype, N_sweep, T_end));
save(filename_out, 'S_star', 'ks_star', 'lam_star', 'band_lo', 'band_hi', 'band_n', 'at_cut', 'conv_cut', ...
    'band_frac', 'ks_cut', 'scheme_run', 'filename_mat', 'FMtype');
fprintf('Saved: %s\n', filename_out);

% =========================================================================
% STEP 5: LATEX TABLE -> tex/
% =========================================================================
if write_table
    if ~exist(tex_folder, 'dir'), mkdir(tex_folder); end
    filename_tex = fullfile(tex_folder, sprintf('SensitivityStar_%s_AllSchemes_N%d_T%d.tex', FMtype, N_sweep, T_end));
    fid = fopen(filename_tex, 'w');
    fprintf(fid, '%% Generated by SensitivityStar.m from %s\n', filename_mat);
    fprintf(fid, '%% Requires booktabs\n');
    fprintf(fid, '\\begin{table}[htbp]\n\\centering\n');
    fprintf(fid, ['\\caption{Sensitivity index at maximal transverse stability, $S_q^{*}$ (dB), of the %s oscillator ', ...
        'per motif and coupling scheme: mean of $20\\log_{10}\\bar{S}_q(\\mathbf{e})$ over the values of $k_s$ at which ', ...
        '$\\lambda^\\perp_{\\max} \\le %.1f\\,\\lambda^\\perp_{\\max}(k_s^{*})$, where $k_s^{*}$ minimizes $\\lambda^\\perp_{\\max}$. ', ...
        '$^{\\dagger}$: $k_s^{*}$ at the end of the computed range, $k_s = %g$. ', ...
        '$^{\\ddagger}$: minimum of $\\lambda^\\perp_{\\max}$ at a value of $k_s$ without a converged index; $k_s^{*}$ is the most stable of the converged values. ', ...
        'Abbreviations as in Table~\\ref{tab:sync_intervals_Rossler}.}\n'], osc_tex{OsID}, band_frac, ks_cut);
    fprintf(fid, '\\label{tab:sens-star-%s}\n', FMtype);
    fprintf(fid, '\\footnotesize\n');
    par_head = strjoin(cellfun(@(s) ['$' s '$'], par_tex_tab, 'UniformOutput', false), ' & ');
    fprintf(fid, '\\begin{tabular}{llccccccc}\n\\toprule\n');
    fprintf(fid, 'Motif & Scheme & $k_s^{*}$ & %s \\\\\n\\midrule\n', par_head);
    rows_t = find(any(has_result(:, scheme_run), 2)).';
    first_ring_row = rows_t(find(rows_t >= first_ring, 1));
    for ti = rows_t
        if ~isempty(first_ring_row) && ti == first_ring_row && any(rows_t < first_ring)
            fprintf(fid, '\\midrule\n');                % Linear | Ring separator
        end
        for CasoBs = scheme_run
            if ~has_result(ti, CasoBs), continue; end
            bc = U.scheme_bc(CasoBs, :);
            cells = strings(1, 6);
            for q = 1:6
                v = S_star(ti, CasoBs, q);
                if isfinite(v), cells(q) = sprintf(['$' num_fmt '$'], v); else, cells(q) = "--"; end
            end
            ks_txt = sprintf(['$' ks_fmt '$'], ks_star(ti, CasoBs));
            if at_cut(ti, CasoBs), ks_txt = [ks_txt '$^{\dagger}$']; end %#ok<AGROW>
            if conv_cut(ti, CasoBs), ks_txt = [ks_txt '$^{\ddagger}$']; end %#ok<AGROW>
            fprintf(fid, '%s & $\\mathbf{b}_%d \\to \\mathbf{c}_%d$ & %s & %s \\\\\n', topo_short{ti}, bc(1), bc(2), ...
                ks_txt, strjoin(cellstr(cells), ' & '));
        end
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
    fclose(fid);
    fprintf('Saved: %s\n', filename_tex);
end

% =========================================================================
% SCRIPT: ShiftSensitivity_xref.m
% DESCRIPTION: One-time conversion of the sensitivity profiles computed
%              with sbar_ij = s_ij * q_j into the normalization
%                  sbar_ij = s_ij * q_j / x_ref
%              of Eq. (eq_normalized_sensitivity), without new runs.
%              x_ref (AttractorScale.m) is one constant per oscillator, so
%              every per-run RMS is divided by x_ref and every index in dB
%              is shifted by -20 log10(x_ref); uncertainties (spreads in
%              dB), convergence and the no-influence test do not change.
%              Files converted (levels 2-4, shared file format):
%                  mat/SensitivityNorm_vs_ks_<Osc>_<topology>_B<C>_N<N>_T<T>.mat
%              dB fields  : norm_dB, norm_dB_all, comp_dB, S_worst, S_mean, S_best
%              linear     : RMS_all (divided by x_ref)
%              new field  : xref
%              A file that already has xref is left as it is, so the script
%              can be run more than once. A copy of every original file is
%              kept in mat/backup_before_xref/.
%              Level 1 (SensitivityIndexMonteCarlo.m) is not converted: run
%              it again (20 runs at one ks). The consolidated file is not
%              converted: run SensitivityNorm_Consolidate.m again.
% =========================================================================
clear; clc;

OsID_run      = [1 2];                 % 1 = Rossler, 2 = Lorenz
mat_folder    = 'mat';
backup_folder = fullfile(mat_folder, 'backup_before_xref');
dB_fields     = {'norm_dB', 'norm_dB_all', 'comp_dB', 'S_worst', 'S_mean', 'S_best'};
lin_fields    = {'RMS_all'};
osc_list      = {'Rossler', 'Lorenz'};

for OsID = OsID_run
    FMtype   = osc_list{OsID};
    xref     = AttractorScale(OsID);
    shift_dB = 20*log10(xref);
    files    = dir(fullfile(mat_folder, sprintf('SensitivityNorm_vs_ks_%s_*_B*_N*_T*.mat', FMtype)));
    files    = files(~contains({files.name}, '_AllTopologies_'));
    fprintf('\n%s: x_ref = %.4f, shift = -%.2f dB, %d profile files\n', FMtype, xref, shift_dB, numel(files));
    n_conv = 0;
    for k = 1:numel(files)
        fn = fullfile(files(k).folder, files(k).name);
        if ismember('xref', who('-file', fn))
            fprintf('  %s: already has xref, left as it is\n', files(k).name);
            continue;
        end
        if ~exist(backup_folder, 'dir'), mkdir(backup_folder); end
        copyfile(fn, fullfile(backup_folder, files(k).name));
        S = load(fn);
        for f = dB_fields
            if isfield(S, f{1}), S.(f{1}) = S.(f{1}) - shift_dB; end      % -Inf and NaN stay
        end
        for f = lin_fields
            if isfield(S, f{1}), S.(f{1}) = S.(f{1}) / xref; end
        end
        S.xref = xref;
        save(fn, '-struct', 'S');
        n_conv = n_conv + 1;
        fprintf('  %s: converted\n', files(k).name);
    end
    fprintf('%s: %d files converted (originals in %s)\n', FMtype, n_conv, backup_folder);
end
fprintf(['\nNext: SensitivityIndexMonteCarlo.m (level 1), SensitivityNorm_vs_ks.m with replot_only = true, ' ...
    'SensitivityNorm_vs_ks_AllSchemes.m (skip_existing = true), SensitivityNorm_Consolidate.m and ' ...
    'SensitivityNorm_PlotAndTable.m for each oscillator.\n']);

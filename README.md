# Topological sensitivity of three-node chaotic network motifs: MATLAB pipeline

Version 2026-11-02. This version includes the attractor-scale normalization:

    sbar_ij = s_ij * q_j / x_ref

This is Eq. (18) of the manuscript. All the scripts are self-contained: the only shared functions are in this folder.

## 1. Requirements

- MATLAB R2020a or later (`tiledlayout`, `exportgraphics`).
- Symbolic Math Toolbox (`syms`, `odeFunction`).
- Parallel Computing Toolbox (optional): used for `parfor` when `use_parallel` / `use_pool` = true.

## 2. Folder layout

Put every `.m` file in one folder (the working folder) and run the scripts from there. The scripts read and write three sub-folders and create them when they are missing:

| Folder | Content |
|---|---|
| `mat/` | Results: LLE, RMS, sensitivity profiles, consolidated files, `AttractorScale_<Osc>.mat` |
| `eps/` | Figures (vector EPS) |
| `tex/` | LaTeX tables |

The symbolic scripts also write generated ODE files in the working folder: `odeMotifCoupled.m`, `odeMotifError.m`, `odeMotifVariational*.m` and `odeMotifSensitivityMismatch*.m`. They are rebuilt automatically, so do not edit them.

Every script has a CONFIGURATION block at the top. Set `OsID` there (1 = Rössler, 2 = Lorenz), together with `topology_idx` / `topology_run` and `CasoBs` / `scheme_run` where they apply. Each stand-alone script starts with `clear; clc; close all`.

## 3. Conventions

- **Oscillators** (`OscillatorType.m`):
  - Rössler: a = b = 0.2, c = 5.7.
  - Lorenz: σ = 10, ρ = 28, β = 2.
- **Topologies** (`topology_idx` 1–13; identical edge lists in every script; each row is `[source target]`):

  | idx | name | edges |
  |---|---|---|
  | 1 | chain | 1→2, 2→3 |
  | 2 | auxiliar | 1→2, 1→3 |
  | 3 | relay | 1↔2, 2↔3 |
  | 4 | mixed_master | 1→2, 2↔3 |
  | 5 | competitive | 1→2, 3→2 |
  | 6 | mixed_slave | 2→1, 2↔3 |
  | 7 | ring_bidir_chain | 1↔2, 2↔3, 3↔1 |
  | 8 | ring_bidir_double_mixed | 1→2, 2↔3, 3↔1 |
  | 9 | ring_bidir_mixed | 1→2, 2↔3, 3→1 |
  | 10 | ring_unidir_chain | 1→2, 2→3, 3→1 |
  | 11 | ring_bidir_auxiliar | 1→2, 1→3, 2↔3 |
  | 12 | ring_unidir_auxiliar | 1→2, 1→3, 2→3 |
  | 13 | ring_bidir_competitive | 1↔2, 2→3, 1→3 |
- **Coupling scheme** `CasoBs = 3(i-1) + j`, which means b_i → c_j (`SISObc.m`: B = e_i, C = e_jᵀ). Coupling: u_k = k_s Σ_{s→k} B C (x_s − x_k).
- **Error system:** e_ij = x_i − x_j. The scripts integrate [x1; e12; e13], and e23 = e13 − e12 is substituted, not integrated.
- **Sensitivity:**
  - Node 1 is nominal; the node-level mismatches are Δq2 and Δq3 (6 parameters: da2 db2 dc2 da3 db3 dc3).
  - The runs start on the synchronization manifold: e(0) = 0, S(0) = 0, x1(0) from an attractor pool, and the node-1 sensitivity rows are held at 0.
  - The index is s_qj(x_i) = 20 log10 rms(sbar_ij), the median over runs.
  - x_ref = vector RMS of the nominal node trajectory over [500, 2500] (`AttractorScale.m`). Rössler ≈ 7.73 (−17.76 dB), Lorenz ≈ 27.67 (−28.84 dB).

## 4. Run order (complete pipeline from scratch)

| Stage | Script | Set | Output |
|---|---|---|---|
| 0 (optional) | `SimulateMotif.m`, `MotifErrorSystem.m` | OsID, topology, CasoBs, ks | Trajectory and error figures (model check) |
| 1 LLE | `MotifLLE_vs_ks_AllTopologies.m` | OsID = 1, then 2 | `mat/LLE_vs_ks_<Osc>_<topology>.mat`, then `LLE_Consolidate.m` runs at the end → `mat/LLE_vs_ks_<Osc>_AllTopologies.mat` |
| 1 LLE | `LLE_PlotAndTable.m` | OsID | `eps/Motif_LLE_vs_ks_AllSchemes_<Osc>_<topology>_N<N>.eps`, `tex/LLE_SyncIntervals_<Osc>_AllTopologies.tex` |
| 2 RMS | `Motif_RMS_vs_ks_AllTopologies.m` | OsID = 1, then 2 | `mat/RMS_vs_ks_<Osc>_<topology>_T2500.mat`, then `RMS_Consolidate.m` runs at the end |
| 2 RMS | `RMS_PlotAndTable.m` | OsID | `eps/Motif_RMS_vs_ks_<Osc>_<topology>_T2500.eps`, `tex/RMS_SyncIntervals_<Osc>_AllTopologies_T2500.tex` |
| 3 Sensitivity, level 1 | `SensitivityIndexMonteCarlo.m` | chain, CasoBs = 1, ks = 2 | Table 2 (`tex/SensitivityVector_*.tex`), `eps/Motif_SensitivityIndex_dB_*`, `_Convergence_*` |
| 3 Sensitivity, level 2 | `SensitivityNorm_vs_ks.m` | chain, CasoBs = 1 and 5 | `eps/SensitivityNorm_vs_ks_<Osc>_chain_B<C>_N25_T2500.eps` |
| 3 Sensitivity, level 3 | `SensitivityNorm_vs_ks_AllSchemes.m` | chain, scheme_run = 1:9 | Table 3 (`tex/TopologicalSensitivity_AllSchemes_*.tex`), `eps/TopologicalSensitivity_AllSchemes_*` |
| 3 Sensitivity, level 4 | `SensitivityNorm_vs_ks_AllTopologies.m` | OsID = 1, then 2 | per-case files, then `SensitivityNorm_Consolidate.m` runs at the end |
| 3 Sensitivity, level 4 | `SensitivityNorm_PlotAndTable.m` | OsID | Table 4 (`tex/`), `eps/TopologicalSensitivity_AllTopologies_<Osc>_b<i>c<j>.eps` |
| 4 Validation | `ValidateSensitivity_vs_ks.m` | chain, CasoBs = 1, da2, r = 0.1 | `eps/ValidationRMS_vs_ks_Rossler_chain_B1_da2_r0.1.eps` |
| 4 Validation | `CompareMotifTopologiesErrors.m` | the two configurations, mismatch | `eps/CompareErrors_*.eps` |

Dependencies between stages:

- Level 2–4 read the LLE file (stage 1, synchronization interval) and the RMS file (stage 2, edge trimming).
- Validation reads the level-2/3 profile.
- Levels 2, 3 and 4 write the same per-case file, `mat/SensitivityNorm_vs_ks_<Osc>_<topology>_B<C>_N<N>_T<T>.mat`. A case computed at one level is reused by the others (`skip_existing = true`).

## 5. Updating existing results to the x_ref normalization (no new runs)

Use this procedure when `mat/` already holds the profiles computed with sbar = s·q_j:

1. `ShiftSensitivity_xref.m`: converts every per-case profile of both oscillators.
   - It calls `AttractorScale.m`.
   - It shifts the dB fields by −20 log10 x_ref, divides `RMS_all` by x_ref and adds the field `xref`.
   - The originals are copied to `mat/backup_before_xref/`. Files that already have `xref` are skipped.
2. `SensitivityIndexMonteCarlo.m`: rerun at ks = 2 (20 runs). This gives Table 2 and its figures.
3. `SensitivityNorm_vs_ks.m` with `replot_only = true`, for CasoBs = 1 and 5: redraws the chain figures without runs.
4. `SensitivityNorm_vs_ks_AllSchemes.m` (`skip_existing = true`): Table 3 and its figure.
5. `SensitivityNorm_Consolidate.m`, then `SensitivityNorm_PlotAndTable.m`, for OsID = 1 and 2: Table 4 and the bar figures.

The validation and comparison figures do not change, so they need no rerun. Levels 3–4 and `SensitivityNorm_Consolidate.m` stop with an error if they meet a profile without `xref`, so old and new results are never mixed.

## 6. Manuscript map

| Manuscript item | File | Script |
|---|---|---|
| Table 1 (Rössler intervals) | `tex/LLE_SyncIntervals_Rossler_AllTopologies.tex`, RMS intervals `tex/RMS_SyncIntervals_Rossler_*` | `LLE_PlotAndTable.m`, `RMS_PlotAndTable.m` |
| Fig. A.1 (LLE charts, 13 panels) | `Motif_LLE_vs_ks_AllSchemes_Rossler_<topology>_N501.eps` | `LLE_PlotAndTable.m` |
| Fig. A.2 (RMS, b3→c3 rings) | `Motif_RMS_vs_ks_Rossler_ring_*_T2500.eps` | `RMS_PlotAndTable.m` |
| Table 2, sensitivity index figure, convergence figure | `SensitivityVector_*.tex`, `Motif_SensitivityIndex_dB_*`, `Motif_SensitivityIndex_Convergence_*` | `SensitivityIndexMonteCarlo.m` |
| Sensitivity norm vs ks (chain b1→c1, b2→c2) | `SensitivityNorm_vs_ks_Rossler_chain_B1/B5_N25_T2500.eps` | `SensitivityNorm_vs_ks.m` |
| Validation figure | `ValidationRMS_vs_ks_Rossler_chain_B1_da2_r0.1.eps` | `ValidateSensitivity_vs_ks.m` |
| Table 3 and chain bar figure | `TopologicalSensitivity_AllSchemes_Rossler_chain_N25_T2500.*` | `SensitivityNorm_vs_ks_AllSchemes.m` |
| Table 4 and Rössler bar figures | `tab:topo-sens-Rossler`, `TopologicalSensitivity_AllTopologies_Rossler_b1c1/b2c2.eps` | `SensitivityNorm_Consolidate.m` + `SensitivityNorm_PlotAndTable.m` |
| Mismatch comparison figure | `CompareErrors_Rossler_chain_CasoBs1_vs_ring_bidir_double_mixed_CasoBs5_dc2_5_10_20pct.eps` | `CompareMotifTopologiesErrors.m` |
| Table B.1 (Lorenz intervals) | `tex/LLE_SyncIntervals_Lorenz_AllTopologies.tex` | `LLE_PlotAndTable.m` (OsID = 2) |
| Lorenz bar figures | `TopologicalSensitivity_AllTopologies_Lorenz_b*c*.eps` | `SensitivityNorm_PlotAndTable.m` (OsID = 2) |

## 7. Numerical settings (defaults of the scripts)

- **LLE** (`MaxTransLyaExp.m`):
  - Method: RK4 + QR on the transverse variational system.
  - Step dt = 0.01 (Rössler) / 0.005 (Lorenz).
  - Averaging time tf = 2000 / 1000, with a node-1 transient of 100 and seed 42.
  - ks grid [0, 5] / [0, 25] with step 0.01.
  - Points within ±0.05 of each zero crossing are recomputed with tf = 2·10⁴ / 10⁴.
- **RMS:**
  - ode45 at 1e-8/1e-10, from the same random initial condition in [0, 1].
  - Window [500, 2500] with Δt = 0.01; synchronized when RMS(e) < 10⁻³.
  - A run stops early when ‖e‖ < 10⁻⁶ for 50 time units.
  - A run counts as diverged when a state exceeds 10⁴.
- **Sensitivity, level 1:**
  - ks = 2, window [200, 1000], dt_sample 0.05, ode45 at 1e-6/1e-8.
  - Median of the runs; stops on convergence (1 dB), at most 100 runs.
- **Sensitivity, levels 2–4:**
  - Window [500, 2500], 25 equally spaced ks over the RMS-trimmed LLE interval.
  - ks with |λ⊥max|·t_tr < 5 are excluded.
  - The same Monte Carlo test as level 1; seeds 42 + r and the same attractor pool at every ks.
  - Lorenz: intervals are cut at ks = 20 (`ks_max_table`), only the first b3→c3 window is used, and linear mixed master b3→c3 is excluded.
- **⟨S̄_q⟩:** mean over the valid ks (internal RMS check passed and convergence met) of 20 log10 sqrt((1/3) Σ_ij S_q(e_ij)²).
- **Validation:**
  - Full network with ode45 at 1e-8/1e-10.
  - Starts on the synchronization manifold from the profile's pool; median of 10 runs.
  - Runs with |x| > 10³ are discarded.

## 8. Files in this package

| File | Role |
|---|---|
| `OscillatorType.m`, `SISObc.m` | Vector fields and parameters; coupling vectors |
| `AttractorScale.m` | x_ref per oscillator (cached in `mat/`) |
| `SimulateMotif.m`, `MotifErrorSystem.m` | Model check: full network and error system |
| `MaxTransLyaExp.m` | LLE engine |
| `MotifLLE_vs_ks_Case.m`, `MotifLLE_vs_ks_AllSchemes.m` | LLE quick looks: one case; one topology with all schemes |
| `MotifLLE_vs_ks_AllTopologies.m`, `LLE_Consolidate.m`, `LLE_PlotAndTable.m` | LLE production run, consolidation, figures and table |
| `Motif_RMS_vs_ks_AllTopologies.m`, `RMS_Consolidate.m`, `RMS_PlotAndTable.m` | RMS production run, consolidation, figures and table |
| `SensitivityIndexMonteCarlo.m` | Sensitivity, level 1 |
| `SensitivityNorm_vs_ks.m` | Sensitivity, level 2 (`replot_only` switch) |
| `SensitivityNorm_vs_ks_AllSchemes.m` | Sensitivity, level 3 |
| `SensitivityNorm_vs_ks_AllTopologies.m`, `SensitivityNorm_Consolidate.m`, `SensitivityNorm_PlotAndTable.m` | Sensitivity, level 4 |
| `ShiftSensitivity_xref.m` | One-time conversion of existing profiles to the x_ref normalization |
| `ValidateSensitivity_vs_ks.m`, `CompareMotifTopologiesErrors.m` | Validation against full nonlinear simulations |

**Removed (superseded or not part of the pipeline):**

- `MotifErrorSystemMismatch.m`
- `VariationalMotif.m` (old copy of `MotifErrorSystem.m` that still integrates e23)
- `MotifSensitivityAnalysis*.m`, `Motifsensitivityanalysis*.m` (early symbolic versions with random initial conditions)
- `SimulateSensitivityMismatch.m`
- `TopologicalSensitivity_Bars.m`, `TopologicalSensitivity_Bars_AllTopologies.m`, `TopologicalSensitivity_CompareMethods.m` (duplicates of level 3/4; `_Bars_AllTopologies` overwrote the figures of `SensitivityNorm_PlotAndTable.m`)
- `Motif_RMS_vs_ks.m`, `Motif_SyncError_RMS_Sweep.m` (early RMS versions with t_transient = 50)
- `S4_Baselines.m`

## 9. Change log

- **2026-11-05:** `PlotValidationRMS.m`: one legend common to the three panels, placed below them (`Layout.Tile = 'south'`, rows of at most 3 entries so the legend fits the figure width) instead of one legend per panel; a panel without predicted influence is labeled with a short text.
- **2026-11-04 (g):** `PlotValidationRMS.m` reads one validation file per mismatch size (`file_base`, `r_list`; missing files are skipped). With several r, the measured RMS(e_ij) of each r is divided by r and drawn against the single first-order prediction x_ref·10^(S_q(e_ij)/20): the curves collapse onto it where the first-order theory holds and fan out where it fails. With one r the figure is the original one. Output `eps/ValidationRMS_vs_ks_<Osc>_<topology>_B<C>_<dq>_r1_5_10_20pct.eps`.
- **2026-11-04 (f):** `CompareMotifTopologiesErrors.m`: any number of configurations in `configs` (rows `[topology, CasoBs, ks]`; default the 2 × 2 design linear chain / ring bidirectional double mixed × b1→c1 / b2→c2 under δc2). `ks = NaN` takes $k_s^*$ from the consolidated file (same rule as `SensitivityStar.m`). Each measured RMS(e_ij) is compared with the first-order prediction from the local index at that $k_s$, RMS ≈ x_ref·10^(S_q(e_ij)/20)·r; console table measured / predicted / ratio, LaTeX table `tex/CompareErrorsTable_<tag>.tex` and `mat/CompareErrors_<tag>.mat`. The B/A difference of ⟨S_q⟩ is no longer printed.
- **2026-11-04 (e):** `PlotValidationRMS.m` and `ValidateSensitivity_vs_ks.m`: with `y_log`, the y-limits span whole decades covering every measured and predicted value, so no curve leaves its panel.
- **2026-11-04 (d):** new `PlotValidationRMS.m`: redraws the validation figure (Fig. 7) from `mat/ValidationRMS_vs_ks_<tag>.mat` without integrating again (same style as `ValidateSensitivity_vs_ks.m`, `y_log` switch).
- **2026-11-04 (c):** `ValidateSensitivity_vs_ks.m`: new switch `y_log` (default true), logarithmic y-axis for Fig. 7.
- **2026-11-04 (b):** `SensitivityNorm_vs_ks.m`: new switch `show_panel_letters` (default false) for the (a)–(f) letters inside the panels of figure 1, which repeated the subfigure letters of Fig. 6. `ValidateSensitivity_vs_ks.m`: `show_band = false` (the band is narrower than the line) and `y_max = []`, which now sets each panel's upper limit from its measured values so no measured curve leaves the panel.
- **2026-11-04:** `SensitivityStar.m` also flags ($^{\ddagger}$, field `conv_cut`) the cases in which the minimum of $\lambda^\perp_{\max}$ falls at a $k_s$ without a converged index; $k_s^*$ is then the most stable of the converged values.
- **2026-11-03 (d):** two new scripts, both run after `SensitivityNorm_Consolidate.m` and need no new runs of the sweeps.
  - `SensitivityStar.m`: sensitivity index at maximal transverse stability, $S_q^*$ = mean of $20\log_{10}\bar{S}_q(\mathbf{e})$ over the band $\lambda^\perp_{\max} \le 0.9\,\min\lambda^\perp_{\max}$ (`band_frac`). Cases whose $k_s^*$ is the last point of an interval cut at the end of the computed range (Rössler 5, Lorenz 20) are flagged. Writes `mat/SensitivityStar_<Osc>_N<N>_T<T>.mat` and `tex/SensitivityStar_<Osc>_AllSchemes_N<N>_T<T>.tex`.
  - `ModalSensitivity_vs_ks.m`: eigenmode prediction of the error sensitivities from master–slave sensitivities at $\kappa = \mu_k k_s$, with the same node-1 pool, window and solver settings as the direct computation, for every diagonalizable motif of one scheme; compares with the direct values. Writes `mat/`, `eps/` and `tex/ModalSensitivity_<Osc>_B<C>_N<N>_T<T>.*`. Run once per scheme (`CasoBs`).
- **2026-11-03 (c):** `SensitivityIndexMonteCarlo.m` colorbar of the component-index figure (Fig. 5): the label is now $s_{q_j}(x_i) = \mathrm{median}_r\,20\log_{10}\mathrm{rms}(\bar{s}_{i,j}^{(r)})$, as in Eq. (eq_median_index); with `stat_predictor = 'SEM'` it reads $20\log_{10}[\mathrm{mean}_r\,\mathrm{rms}(\bar{s}_{i,j}^{(r)})]$, which is what `indexStat` computes.
- **2026-11-03 (b):** `SensitivityIndexMonteCarlo.m` convergence figure: title uses $N_M$ and "realizations" (lowercase statistic), x-axis "Number of realizations $n$", matching the manuscript.
- **2026-11-03:** Scheme labels in figure titles, axis labels and generated LaTeX tables are now bold, $\mathbf{b}_p \to \mathbf{c}_q$, matching the manuscript (10 scripts, display strings only).
- **2026-11-02 (d):** `SensitivityNorm_PlotAndTable.m` writes the definition of the estimator gap into the table caption.
- **2026-11-02 (c):**
  - Topology 8 is now called "Ring bidirectional double mixed" in every script, matching its key `ring_bidir_double_mixed`, the tables and Fig. 1.
  - Bar figures (`SensitivityNorm_vs_ks_AllSchemes.m`, `SensitivityNorm_PlotAndTable.m`): the x-axis comes from the new local function `barLimits`.
    - The zero line is included and the ends are multiples of 5 dB.
    - A margin of 4 % of the range is left on each side, or 14 % when the values are written, so no bar or label touches the frame.
  - The all-schemes figure is 120 px tall per scheme (was 70), so the value labels of neighbouring bars no longer overlap.
- **2026-11-02 (b):**
  - `topoSensIndex` in `SensitivityNorm_Consolidate.m` and `SensitivityNorm_vs_ks_AllSchemes.m`: a vector that is structurally zero at a valid ks (NaN norm while other vectors are finite) now counts as 0 instead of making ⟨S̄_q⟩ NaN. This fixes Lorenz chain b1→c1 δa2, which was shown as "--" (exact first-order result: s(e23) = 0, see section 10).
  - A ks where every vector is NaN is left out.
  - `SensitivityNorm_vs_ks.m`: the console line of `replot_only` no longer says the mat file was saved.
- **2026-11-02:**
  - Attractor-scale normalization: new `AttractorScale.m` and `ShiftSensitivity_xref.m`.
  - The four level scripts now divide by x_ref and save `xref`.
  - Levels 3–4 and `SensitivityNorm_Consolidate.m` check for `xref`.
  - `ValidateSensitivity_vs_ks.m` multiplies the prediction by x_ref.
  - New `replot_only` switch in `SensitivityNorm_vs_ks.m`.
- **2026-11-02:** `RMS_PlotAndTable.m` uses a log-scale RMS axis with the sync_tol line and shaded diverged runs.
- **2026-11-01:**
  - `MotifErrorSystem.m`: e23 = e13 − e12 substituted in every row, [x1; e12; e13] integrated, and e23 rebuilt for the output.
  - `LLE_Consolidate.m`: saves `LLE_coarse_array`, `refined_mask_array` and `opts_refine_all`.
  - `SimulateMotif.m`: blow-up value.
  - `MaxTransLyaExp.m` and `MotifLLE_vs_ks_Case.m`: header and label fixes.

## 10. Known points

- Lorenz, b1→c1 (x coupling), mismatch δσ at node 2: to first order the node-2 sensitivity is ẋ1/(σ k_s), a time shift of the trajectory. Node 3 receives the same signal, so s(e12) = s(e13) and s(e23) = 0 exactly; the remaining e23 is second order. The full-network simulation confirms it.
- The estimator gap (mean of the per-ks medians vs median of the per-run indices) reaches 3 dB for Lorenz, mostly under b3→c3, because the per-run distributions are heavy tailed. The Rössler gap stays below 0.22 dB.
- The divergence bound is 10³ in the sensitivity and validation scripts and 10⁴ in the RMS scripts.

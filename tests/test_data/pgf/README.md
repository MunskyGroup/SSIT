# PGF reference fixtures (phase 0)

Reference data for the native MATLAB PGF solver (`src/+ssit/+pgf`). Generated from Edward's Julia
implementation by `PGFfromEdward/export_pgf_fixtures.jl`:

```
cd PGFfromEdward
julia --project=. -t 4 --startup-file=no export_pgf_fixtures.jl      # writes into this folder
```

All `.mat` files are v7.3 (HDF5) and load with MATLAB's `load`. Integer data (multi-indices, indices) are
stored as doubles. Strings are cell arrays. Indices are 1-based.

## Conventions

- `grid` is N x d, one node per row. `domain` is 2 x d: row 1 holds the upper bounds (always 1), row 2 the lower bounds.
- Nodes are Chebyshev–Gauss–Lobatto points (Chebyshev extrema) from SmolyakApprox.jl.
  - Level i has m(i) = 1 (i = 1) or 2^(i-1)+1 nodes. Each level contributes only the nodes that are new relative to
    the previous level.
  - Node ordering follows `smolyak_grid`: blocks in multi-index order, and inside a block the first dimension varies
    fastest.
- An integer `mu` uses the isotropic multi-index generator and a vector `mu` the anisotropic one. Both appear in
  `pgf_grids.mat`. For equal entries they give the same multi-index, row order included (checked for d ≤ 4, μ ≤ 3).
- `D(:,:,i)` is d/dz_i at the nodes. `H(:,:,a,b)` is d²/(dz_a dz_b). Both come from the Smolyak Chebyshev
  interpolant, built from the basis matrix times `inv_interp_mat`.
- `one_index` is the row of `grid` equal to all ones. It is 0 when no such node exists.

## Files

| File | Contents | Reference solution |
|---|---|---|
| `pgf_grids.mat` | 5 grid cases (d = 1–4, isotropic and anisotropic): grid, multi-index, `inv_interp_mat`, `D`, `H` (d ≤ 3) | Edward's code (exact port target, about 1e-12) |
| `pgf_birthdeath.mat` | 0→X (k), X→0 (g·x1); L, G(t) | `expm`. Matches the analytic Poisson to 6e-14 |
| `pgf_chain2.mat` | 0→x1, x1→x2, decays (Edward's `rhs_chain_pgf!`); L, G(t) | `expm`. Analytic product-Poisson given for information; discretisation error 1.5e-5 at μ=[3,3] |
| `pgf_annihilation.mat` | 0→X (k), 2X→0 with propensity c·x1·(x1−1); domain [−1, 1] | `expm`. Matches FSP (`p_fsp`, `Gfsp`) to 3e-13 |
| `pgf_grdusp1.mat` | GR→DUSP1 model: genes C/D enumerated as blocks [g0; g1], PGF over A, B, E, F; μ=[4,2,2,2] on [0.95,1]^4. For `truth` and `ones`: `L_pre*v_test`, `L_post*v_test`, steady state, 9 snapshots, MSE objective. Also the empirical data PGF (`ref_pgf`, `ref_pgf_sq`) | Rodas5P at 1e-12 (see below). `gtot_qndf` = Edward's production settings |
| `grdusp1_steady_counts.csv` | The 5000-cell SSA data behind `ref_pgf`, in SSIT's CSV layout (`time,A,...,F`, hours) | — |
| `bd_steady_counts.csv` | Phase 4: 2000 draws from Poisson(10) at t = 0 (steady state of 0→X at k = 10, X→0 at g = 1); `time,x1` | Exact MLE of k = sample mean |
| `bd_timecourse_counts.csv` | Phase 4: 500 draws each at t = 1, 2, 4, 8 from Poisson(k/g·(1 − e^(−gt))), k = 10, g = 0.5, x(0) = 0 | Exact Poisson MLE (computed in `pgfFitTest`) |

Objective check: unweighted MSE at truth = 1.0398e-5 and at all-ones = 0.22440. Edward's QNDF probe gives
1.043e-5 and 0.2244.

## Findings that constrain the MATLAB design

1. **The all-ones node exists only when max(μ) ≥ d.** Example: `d3_aniso212` has no z = 1 node. The method relies
   on it for normalisation (G(1) = 1, used in the steady-state solve) and for moments. Grid options must enforce
   max(μ) ≥ d, or the code must evaluate z = 1 through the interpolant.
2. **Second-order terms in a PGF species need the full domain [−1, 1] for that species.** For 2X→0 the PDE term
   is c(1−z²)G''. On a truncated [zmin, 1] the collocation operator has eigenvalues with large positive real part:
   at μ=4, +19 on [0.8, 1] and +114 on [0.95, 1], and these grow with refinement. On [−1, 1] it is stable.
   First-order operators are stable on any [zmin, 1]. So a propensity of degree ≥ 2 in a PGF species is allowed
   only if that species' domain is [−1, 1]; otherwise the species must be enumerated (as Edward does with B+C,
   where C is enumerated) or the model rejected.
3. **The operators can be strongly non-normal.**
   - At all-ones parameters, `L_post` is spectrally stable (max Re λ ≈ 0). But a relative perturbation of 1e-12
     moves its eigenvalues to Re λ = +0.9, and the computed ‖expm(L_post·25 h)‖ reaches about 1e42.
   - Do **not** use `expm`/`expmv` for the solves; use an implicit stiff integrator (ode15s/ode23s).
   - Stiff solvers stay bounded and agree, but only to about 1e-7 at all-ones, against about 1e-11 at truth.
     Test tolerances for the MATLAB port should reflect this.

## Phase 4 (fitting) notes

- The birth–death CSVs come from numpy (`default_rng(20261109)`), written by the session's mirror script
  `verify_phase4_bd.py`. That script also predicts the PGF estimates the tests check.
- `gls` loss on the GR data: with the full covariance of the empirical PGF, χ² at the true parameters is calibrated
  (311 against 255 effective degrees of freedom) only when the PGF solution is accurate. That holds at
  μ = [5,3,3,3], but not at Edward's [4,2,2,2], where the grid error (3.75e-4 RMS) dominates the small eigendirections.
  This is why `errorFloor = 'auto'` adds the RMS change under μ + 1 to the noise model.

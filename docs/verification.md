# Verification evidence

Every number below is a real measurement from the tested runtime
**GNU Octave 11.3.0** (Windows x86_64), produced by
`tests/test_heatwave.m` (9/9 passed, 23.3 s) and `run_results.m`.
Machine-readable copies: `docs/results/run_summary.json`.

## Test 1 — defaults (`test_defaults`)

`heatwave()` with all defaults: `nx=50`, `nSteps=556`, `r=0.4496`
(auto-selected below the stability limit), lands exactly on
`tFinal = 0.1`, stores every step, initial row equals `sin(pi x)` to
roundoff. Passed.

## Test 2 — analytic solution self-checks (`test_exact_solution`)

- `u(x,0) = sin(pi x)` to `1e-14` (single mode) and to `1e-14` for a
  3-mode series `c = [0.5, -0.3, 0.2]`.
- Independent PDE check: centered finite-difference `du/dt` of
  `exact_solution` at `(x,t) = (0.5, 0.4)` vs the analytic rate
  `-pi^2 u`: residual **3.09e-08** (consistent with the `O(ht^2)`
  finite difference, `ht = 1e-4`).
- Per-mode decay ratios match `exp(-k^2 pi^2 t)` for `k = 1, 3` to
  `1e-14`.
- Superposed steady profile (`uLeft=1, uRight=3`) keeps `u(0,t)=1`,
  `u(1,t)=3` and equals `1+2x+0.7 sin(pi x) e^{-pi^2 t}` to `1e-12` at
  `t = 0, 0.2, 1.7`. Passed.

## Test 3 — FTCS vs analytic decay (`test_ftcs_vs_analytic`)

`nx=40` (`dx=0.025`), `nSteps=720` (`dt=2.778e-4`, `r=0.444`),
`u0 = sin(pi x)`, zero Dirichlet ends. Max-norm error against
`sin(pi x) e^{-pi^2 t}`:

| t | max abs error |
|---|---|
| 0.05 | 2.584e-04 |
| 0.10 | 3.155e-04 |
| 0.20 | 2.350e-04 |

Tolerance asserted: `< 6e-4`, which is consistent with the
`O(dt) + O(dx^2)` truncation error at this resolution (measured values
are roughly half the tolerance). Passed.

## Test 4 — implicit at a time step where FTCS cannot run (`test_implicit_large_dt`)

`nx=40`, `dt=0.01`, `tFinal=0.2`: mesh ratio `r = 16`, i.e. **32x** the
FTCS stability limit `dx^2/(2 alpha) = 3.125e-4`.

- Backward Euler + Thomas: max error vs analytic **1.344e-02**
  (asserted `< 2.5e-2`); this is honest first-order accuracy at a
  deliberately coarse step — the point is stability, not extra accuracy.
- FTCS at the same `dt` with `'stabilityGuard',false` and initial data
  `sin(pi x) + 0.05 sin(60 pi x)`: final max norm **3.36e+28**
  (diverged, asserted `> 1e6`). A run with the pure `sin(pi x)` initial
  condition also diverges (measured `2.8e+19`) because roundoff alone
  seeds the unstable high modes; the perturbed case is used in the test
  because its growth is deterministic. Passed.

## Test 5 — stability guard (`test_stability_guard`)

`nx=10` (`dx=0.1`, limit `dx^2/(2 alpha) = 0.005`):

- `dt=0.006` raises `heatwave:unstableStep` (guard fired).
- `dt=0.005` (exactly at the limit) is accepted and used **unclipped**
  (`sol.info.dt == 0.005` to machine precision). Passed.

## Test 6 — spatial convergence (`test_convergence_dx`)

Fixed `r` near 0.45 (auto step count), `tFinal=0.1`, analytic reference:

| nx | dx | nSteps | max error |
|---|---|---|---|
| 10 | 1.000e-01 | 23 | 4.943e-03 |
| 20 | 5.000e-02 | 89 | 1.288e-03 |
| 40 | 2.500e-02 | 356 | 3.211e-04 |
| 80 | 1.250e-02 | 1423 | 8.031e-05 |
| 160 | 6.250e-03 | 5689 | 2.009e-05 |

Measured order (log-log fit): **1.989** (asserted 1.8–2.2). With `r`
fixed, both the `O(dt)` and `O(dx^2)` components scale like `dx^2`, so
this is the composite second-order behavior of FTCS. Passed.

## Test 7 — temporal convergence (`test_convergence_dt`)

Fixed `nx=200`, `tFinal=0.05`, `dt` halved from `r=0.5` to `r=1/32`.
Measured by self-convergence `max|U(dt) - U(dt/2)|`:

| dt | \|U(dt)-U(dt/2)\| | \|U-exact\| |
|---|---|---|
| 1.250e-05 | 9.293e-06 | 1.239e-05 |
| 6.250e-06 | 4.646e-06 | 3.097e-06 |
| 3.125e-06 | 2.323e-06 | 1.549e-06 |
| 1.563e-06 | 1.161e-06 | 3.872e-06 |

Measured order: **1.000** (asserted 0.8–1.2). The self-convergence
column halves exactly with `dt`. The raw analytic-error column is
deliberately reported as well: it is *non-monotone* (minimum near
`r = 1/6`, rising again toward the spatial floor `5.6e-6` at `r = 1/32`)
because FTCS's leading time and space error terms partially cancel —
see `docs/model.md` section 3. Self-convergence is insensitive to that
cancellation and is the honest measure of the temporal order. Passed.

## Test 8 — insulated-slab mean conservation (`test_neumann_mean`)

`nx=60`, `nSteps=4000`, `tFinal=0.5`, `bc='neumann'`,
`u0 = 1 + 0.5 sin(pi x)`:

| scheme | max drift of node mean | final spread about mean |
|---|---|---|
| ftcs | **1.332e-15** | 1.077e-09 |
| implicit | **1.787e-13** | 1.180e-09 |

The mean is conserved to roundoff (asserted `< 1e-10`), matching the
exact telescoping property of the ghost-cell update
(`docs/model.md` section 4). The final profile is flat to ~1e-9 because
all cos-modes present in the initial condition decay at least as fast as
`e^{-4 pi^2 t}`. Passed.

## Test 9 — input validation (`test_validation`)

23 malformed-input cases (bad `nx`/`L`/`alpha`/`tFinal`/`dt`, unknown or
duplicate options, odd argument count, wrong-length or NaN `u0`,
conflicting `dt`+`nSteps`, non-dividing `dt`, bad `scheme`/`bc`, bad
`captureEvery`, non-logical guard, and 5 malformed `exact_solution`
inputs) all raise `heatwave:*` / `exact_solution:*` errors. A valid
vector `u0` is stored bit-exact. Passed.

## Results artifacts (`run_results.m`, same runtime)

- `evolution.png` — FTCS profiles, `nx=120`, `nSteps=16000`, `r=0.450`,
  13 stored times, pure white figure interior, no grid lines.
- `evolution_surface.png` — mesh of `u(x,t)` over the stored states.
- `convergence_dx.png`, `convergence_dt.png` — the tables above with
  slope-2 / slope-1 reference lines.
- `error_final.png` — pointwise `|FTCS - analytic|` at `t=0.5` for
  `nx=120`: max **3.45e-06**.
- `frames/frame_001..062.png` — one PNG per stored time level of that
  run; `animation.gif` (62 frames, 80 ms each) is assembled from these
  exact files by `tools/make_gif.py` (Python + Pillow). No frame is
  drawn from anything but a computed state.
- `run_summary.json` — the numbers above in machine-readable form.

## Limitations

- One spatial dimension; uniform grid; node-centered unknowns.
- Constant `alpha` (no temperature-dependent diffusivity).
- Linear problem: no convection, sources, or radiation.
- Boundary values are constant in time.
- The Thomas solve omits pivoting — valid here because the backward
  Euler matrices are strictly diagonally dominant by construction; it is
  not a general-purpose tridiagonal solver.
- MATLAB compatibility is a design goal; the only tested runtime is
  GNU Octave 11.3.0.

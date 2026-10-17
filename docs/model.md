# Model, discretizations, and analytic solution

This document states exactly what Heatwave solves, how it discretizes it,
why the time-step limit exists, and what the analytic reference is. All
measurements quoted here were produced on the tested runtime,
**GNU Octave 11.3.0** (Windows, x86_64).

## 1. Governing equation

Transient heat (mass) diffusion in a one-dimensional slab of length `L`:

```
du/dt = alpha * d2u/dx2,     0 <= x <= L,  0 <= t <= tFinal
```

with a prescribed initial temperature `u(x,0) = u0(x)` and either:

- **Dirichlet ends**: `u(0,t) = uLeft`, `u(L,t) = uRight` (fixed end
  temperatures), or
- **Neumann (insulated) ends**: `du/dx = 0` at `x=0` and `x=L`, so the slab
  exchanges no heat with its surroundings and the mean temperature is
  conserved.

### Symbols and units

| Symbol | Option / field | Meaning (illustrative SI units) |
|---|---|---|
| `u` | state | temperature (K) or concentration (kg/m^3) |
| `x` | `sol.x` | position along the slab (m) |
| `t` | `sol.t` | time (s) |
| `alpha` | `'alpha'` | thermal diffusivity (m^2/s) |
| `L` | `'L'` | slab length (m) |
| `dx` | `sol.info.dx` | grid spacing `L/nx` (m) |
| `dt` | `sol.info.dt` | time step (s) |
| `r` | `sol.info.r` | mesh ratio `alpha*dt/dx^2` (-) |

## 2. Grid

Uniform node-centered grid `x_i = i*dx`, `i = 0..nx`, `dx = L/nx`
(`nx+1` nodes). All dependent variables are stored at the nodes.

## 3. Explicit FTCS scheme

Forward Euler in time, centered second difference in space. For interior
nodes:

```
u_i^{n+1} = u_i^n + r * (u_{i+1}^n - 2 u_i^n + u_{i-1}^n),
r = alpha*dt/dx^2
```

Boundary rows:

- Dirichlet: `u_0^{n+1} = uLeft`, `u_nx^{n+1} = uRight` (imposed, not
  updated by the PDE).
- Neumann, ghost cells: with fictitious nodes `u_{-1} = u_0`,
  `u_{nx+1} = u_{nx}` (zero gradient), the end nodes become
  `u_0^{n+1} = u_0^n + r (u_1^n - u_0^n)` and symmetrically at `x = L`.

Local truncation error: `O(dt) + O(dx^2)`. The implementation is fully
vectorized: one slice expression updates all interior nodes per step.

### Stability condition (derivation sketch)

Apply von Neumann analysis to the interior update: a Fourier mode
`u_i^n = G^n e^{j i theta}` obeys

```
G = 1 - 4 r sin^2(theta/2)
```

The scheme is stable iff `|G| <= 1` for every resolvable `theta`. Since
`0 <= sin^2(theta/2) <= 1`, this requires `|1 - 4r| <= 1`, i.e.

```
r <= 1/2    <=>    dt <= dx^2 / (2 alpha).
```

Heatwave enforces this exactly: for `'scheme','ftcs'` with the default
`'stabilityGuard',true`, a time step above the limit raises a
`heatwave:unstableStep` error. The step is **never silently clipped**.
The guard can only be switched off by passing `'stabilityGuard',false`
explicitly, which the test suite uses to demonstrate the instability.

### Measured note: partial error cancellation near r = 1/6

For a smooth fundamental mode the leading `O(dt)` and `O(dx^2)` error
terms of FTCS have opposite signs. Expanding the per-step amplification
`G = 1 - 4r sin^2(theta/2)` against `exp(-theta^2 r)` for
`theta = pi*dx/L` gives an excess decay rate proportional to
`(1/12 - r/2)`, which vanishes at `r = 1/6`. Measured on this runtime at
`nx = 200`, `tFinal = 0.05`, the analytic error `|U - u_exact|` is
non-monotone in `dt`: `1.24e-5, 3.10e-6, 1.55e-6, 3.87e-6, 5.61e-6` for
`r = 0.5, 0.25, 0.125, 0.0625, 0.016` (minimum near `r = 1/6`, then it
rises back to the spatial floor). This is exactly why the temporal
convergence test measures order by **self-convergence**
(`|U(dt) - U(dt/2)|`), which is unaffected by the cancellation, instead
of by raw analytic error; see `docs/verification.md`.

## 4. Implicit scheme (backward Euler + own Thomas solver)

Backward Euler in time, centered in space:

```
-r u_{i-1}^{n+1} + (1 + 2r) u_i^{n+1} - r u_{i+1}^{n+1} = u_i^n
```

- Dirichlet: the endpoint values are known, so the system holds for the
  `nx-1` interior unknowns with diagonal `1+2r` and off-diagonals `-r`;
  the boundary coefficients `r*uLeft`, `r*uRight` move to the RHS.
- Neumann: ghost-cell rows reduce the end equations to
  `(1+r) u_0^{n+1} - r u_1^{n+1} = u_0^n` (and symmetrically), giving a
  `nx+1` system.

The tridiagonal system is solved by a **hand-written Thomas algorithm**
(forward elimination without pivoting, then back substitution) in
`src/heatwave.m`. No pivoting is needed because both matrix variants are
strictly diagonally dominant (`1+2r > 2r`; `1+r > r`), which also
guarantees the elimination cannot divide by zero.

Backward Euler is unconditionally stable: its amplification factor
`1/(1 + 4r sin^2(theta/2))` satisfies `|G| <= 1` for every `r >= 0`,
so large steps degrade accuracy (first order in `dt`) but never blow up.
Local truncation error: `O(dt) + O(dx^2)`.

### Discrete mean conservation (Neumann case)

Let `A` be the ghost-cell Neumann Laplacian matrix (rows `[-1 1 ...]`,
`[1 -2 1]`, ..., `[..., 1 -1]`). Every column of `A` sums to zero, so for
both schemes `(I + rA) u^{n+1} = u^n` or the FTCS equivalent implies

```
sum_i u_i^{n+1} = sum_i u_i^n      (exactly, in exact arithmetic)
```

i.e. the arithmetic node mean is conserved; on the tested runtime the
drift over 4000 steps is at roundoff level: `1.3e-15` (FTCS) and
`1.8e-13` (implicit).

## 5. Analytic solution

Separation of variables on `[0, L]` with Dirichlet ends
`u(0,t)=uLeft`, `u(L,t)=uRight`: the steady state is the linear profile
`S(x) = uLeft + (uRight - uLeft) x/L`, and the transient part is a
decaying sine series. For a single mode of amplitude `A` and mode number
`m` (the documented reference case uses `A = 1`, `m = 1`, zero ends):

```
u(x,t) = S(x) + A sin(m pi x / L) exp(-alpha (m pi / L)^2 t)
```

and for a general initial condition (finite sine series with
coefficients `c_k`):

```
u(x,t) = S(x) + sum_k c_k sin(k pi x / L) exp(-alpha (k pi / L)^2 t)
```

Each mode decays exponentially with rate `alpha (k pi / L)^2`. For
insulated ends the eigenfunctions are cosines and the `k = 0` mode is
constant, which is the analytic counterpart of discrete mean
conservation.

`src/exact_solution.m` implements the single-mode and series forms with
the linear steady profile; the suite verifies it independently by a
finite-difference PDE residual check (measured `3.1e-8` at `(x,t) =
(0.5, 0.4)`), a series identity at `t = 0`, per-mode decay ratios, and
boundary values.

## 6. MATLAB / GNU Octave compatibility statement

- The `src/` and `tests/` code uses only core-language constructs
  (numeric arrays, function handles, `error`/`assert`, `struct`); no
  toolbox, no graphics in tests, no platform-specific calls.
- **Tested runtime: GNU Octave 11.3.0** on Windows x86_64 (locally) and
  the Octave release shipped with `ubuntu-latest` in GitHub Actions.
- MATLAB compatibility is a design goal, not a tested claim: the sources
  are written against the common MATLAB/Octave core language
  (R2016b-and-later conventions), but no MATLAB license was available to
  this project for verification.

## 7. References (by title)

- *Numerical Solution of Partial Differential Equations: An Introduction*
  (K. W. Morton, D. F. Mayers) — FTCS stability and accuracy.
- *Finite Difference Schemes and Partial Differential Equations*
  (J. C. Strikwerda) — von Neumann analysis.
- *Numerical Recipes: The Art of Scientific Computing* (W. H. Press et
  al.) — tridiagonal (Thomas) solve.
- *Heat Conduction* (D. W. Hahn, M. N. Ozisik) — separation-of-variables
  series solutions.

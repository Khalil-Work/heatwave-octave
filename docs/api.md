# API reference

All numerical entry points live in `src/`. Tested runtime:
**GNU Octave 11.3.0**; MATLAB compatibility is a stated goal (core
language only), see `docs/model.md` section 6.

## `heatwave` — transient diffusion solver

```matlab
sol = heatwave()                                  % documented defaults
sol = heatwave('nx', 100, 'tFinal', 0.2, ...)     % name/value options
```

### Options

| Name | Default | Meaning |
|---|---|---|
| `'L'` | `1` | slab length (positive scalar) |
| `'nx'` | `50` | number of grid intervals; grid has `nx+1` nodes (integer >= 2) |
| `'alpha'` | `1` | diffusivity (positive scalar) |
| `'tFinal'` | `0.1` | final time (positive scalar) |
| `'dt'` | *(auto)* | time step; `tFinal/dt` must be an integer within 1e-9 relative, otherwise an error suggests `'nSteps'` |
| `'nSteps'` | *(auto)* | step count; `dt = tFinal/nSteps`. Mutually exclusive with `'dt'` |
| *(neither given)* | — | FTCS: smallest count with `r <= 0.45`; implicit: 200 steps |
| `'scheme'` | `'ftcs'` | `'ftcs'` (explicit) or `'implicit'` (backward Euler + Thomas) |
| `'bc'` | `'dirichlet'` | `'dirichlet'` (`'uLeft'`/`'uRight'`) or `'neumann'` (insulated, ghost cells) |
| `'u0'` | *(auto)* | function handle `u0(x)` or vector of `nx+1` values; empty selects `sin(pi x/L)` |
| `'uLeft'` | `0` | Dirichlet value at `x=0` (ignored for Neumann) |
| `'uRight'` | `0` | Dirichlet value at `x=L` (ignored for Neumann) |
| `'captureEvery'` | `1` | store every N-th step plus step 0 and the final step |
| `'stabilityGuard'` | `true` | if true, FTCS **errors** (never clips) when `dt > dx^2/(2 alpha)` |

### Returns (`sol`)

| Field | Size | Content |
|---|---|---|
| `sol.x` | `nx+1 x 1` | node coordinates |
| `sol.t` | `nTimes x 1` | stored times (`0 ... tFinal`) |
| `sol.U` | `nTimes x nx+1` | stored states; row `k` is `u` at `sol.t(k)` |
| `sol.info` | struct | all parameters above plus `dx`, `r`, `dtMaxStable = dx^2/(2 alpha)`, `guardApplied` |

### Error identifiers

`heatwave:usage`, `heatwave:invalidOption`,
`heatwave:invalidInput`, `heatwave:unstableStep`.

## `exact_solution` — analytic reference

```matlab
u = exact_solution(x, t)                       % sin(pi x/L) exp(-alpha pi^2 t/L^2)
u = exact_solution(x, t, 'L', L, 'alpha', alpha, ...)
```

| Name | Default | Meaning |
|---|---|---|
| `'L'` | `1` | slab length |
| `'alpha'` | `1` | diffusivity |
| `'mode'` | `1` | sine mode number `m` (positive integer) |
| `'amplitude'` | `1` | amplitude `A` |
| `'coeffs'` | `[]` | if non-empty, overrides mode/amplitude: `u = S(x) + sum_k c(k) sin(k pi x/L) exp(-alpha (k pi/L)^2 t)` |
| `'uLeft'` | `0` | steady-state value at `x=0`; `S(x)` is the linear steady profile |
| `'uRight'` | `0` | steady-state value at `x=L` |

`x` is a vector (orientation preserved), `t` a non-negative scalar.
Errors: `exact_solution:usage`, `exact_solution:invalidOption`,
`exact_solution:invalidInput`.

## `tests/test_heatwave.m` — verification suite

```bash
octave --no-gui --eval "run('tests/test_heatwave.m')"   # from repo root
```

Function file; prints per-test measured numbers, exits `0` on full pass
and `1` on any failure. No graphics, no toolboxes (CI-safe).

## `run_results.m` — regenerate all figures and numbers

```bash
octave --no-gui --eval "run('run_results.m')"    # writes docs/results/*
python tools/make_gif.py                          # assembles animation.gif
```

Overwrites `docs/results/` figures, `run_summary.json`, and
`frames/frame_*.png` (stale frames are removed first). On the tested
Windows install it uses the gnuplot graphics toolkit with
`'visible','off'` figures and `print -dpng` (offscreen Qt PNG export
fails there); that mechanism is part of the documented behavior.

## `tools/make_gif.py` — GIF assembly

```bash
python tools/make_gif.py [--frames DIR] [--out FILE] [--duration MS]
```

Sorts `frame_*.png`, quantizes each to an adaptive 256-color palette,
and writes a looping GIF (default 80 ms/frame). Requires Pillow. Used
because Octave-side animated GIF export is unreliable across builds;
see the provenance caption in the README.

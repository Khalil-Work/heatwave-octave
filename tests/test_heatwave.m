function test_heatwave()
%TEST_HEATWAVE Verification suite for the Heatwave solvers (GNU Octave).
%
%   Exit status: 0 when every test passes, 1 otherwise.  Each test prints
%   its measured numbers, which are reproduced in docs/verification.md.
%
%   Run from the repository root with either:
%     octave --no-gui --eval "run('tests/test_heatwave.m')"
%   or directly:
%     octave --no-gui tests/test_heatwave.m
%
%   Tested runtime: GNU Octave 11.3.0 (Windows x86_64).  The suite uses
%   only core-language features (no toolboxes, no graphics) so it also
%   runs under MATLAB and on headless CI.

more off;

% --------------------------------------------------------------- locate src
here = fileparts(mfilename('fullpath'));
src_candidates = {};
if ~isempty(here)
  src_candidates{end+1} = fullfile(fileparts(here), 'src'); %#ok<SAGROW>
end
src_candidates{end+1} = fullfile(pwd, 'src'); %#ok<SAGROW>
srcdir = '';
for k = 1:numel(src_candidates)
  if exist(src_candidates{k}, 'dir')
    srcdir = src_candidates{k};
    break;
  end
end
if isempty(srcdir)
  error('test_heatwave:setup', ...
        'cannot locate the src/ directory; run from the repository root.');
end
addpath(srcdir);
fprintf('Heatwave verification suite (GNU Octave %s)\n', version());
fprintf('src: %s\n\n', srcdir);

% ------------------------------------------------------------------ driver
tests = { ...
  @test_defaults, ...
  @test_exact_solution, ...
  @test_ftcs_vs_analytic, ...
  @test_implicit_large_dt, ...
  @test_stability_guard, ...
  @test_convergence_dx, ...
  @test_convergence_dt, ...
  @test_neumann_mean, ...
  @test_validation};

n_fail = 0;
t_start = tic;
for k = 1:numel(tests)
  name = func2str(tests{k});
  try
    feval(tests{k});
    fprintf('PASS  %s\n\n', name);
  catch err
    n_fail = n_fail + 1;
    fprintf(2, 'FAIL  %s\n       id: %s\n       msg: %s\n\n', ...
            name, err.identifier, err.message);
  end
end
elapsed = toc(t_start);
fprintf('suite: %d/%d passed in %.1f s\n', numel(tests)-n_fail, ...
        numel(tests), elapsed);

if n_fail == 0
  fprintf('ALL TESTS PASSED\n');
  exit(0);
else
  fprintf(2, '%d TEST(S) FAILED\n', n_fail);
  exit(1);
end
end

% ============================================================================
function test_defaults()
% Defaults run, land exactly on tFinal, and respect the FTCS limit.
  s = heatwave();
  assert(s.info.r <= 0.5 + 1e-12);
  assert(abs(s.t(end) - 0.1) < 1e-12);
  assert(isequal(size(s.U), [s.info.nSteps + 1, s.info.nx + 1]));
  assert(max(abs(s.U(1,:) - sin(pi*s.x).')) < 1e-15);
  assert(max(abs(s.U(end,:) - exact_solution(s.x, 0.1).')) < 1e-2);
  fprintf('    defaults: nx=%d nSteps=%d r=%.4f dt=%.3e\n', ...
          s.info.nx, s.info.nSteps, s.info.r, s.info.dt);
end

% ============================================================================
function test_exact_solution()
% Analytic formula self-checks (independent finite-difference residual,
% series identity, boundary values, steady-state superposition).
  % documented default case: u = sin(pi x) exp(-pi^2 t) on [0,1]
  x = linspace(0, 1, 201);
  assert(max(abs(exact_solution(x, 0) - sin(pi*x))) < 1e-14);
  % PDE residual at (x,t) = (0.5, 0.4): u_t - alpha u_xx ~ 0 via FD in t
  % (u_xx = -pi^2 u for the eigenmode, so the check uses that identity)
  t0 = 0.4; ht = 1e-4;
  ut_fd = (exact_solution(0.5, t0+ht) - exact_solution(0.5, t0-ht)) / (2*ht);
  ut_true = -pi^2 * exact_solution(0.5, t0);
  assert(abs(ut_fd - ut_true) < 1e-6);
  fprintf('    PDE residual at (0.5,0.4): %.3e\n', abs(ut_fd - ut_true));
  % general mode and amplitude
  assert(max(abs(exact_solution(x, 0, 'mode', 2, 'amplitude', 0.5) ...
              - 0.5*sin(2*pi*x))) < 1e-14);
  % decaying sine series: u(x,0) reproduces the superposed modes
  c = [0.5, -0.3, 0.2];
  u0_series = c(1)*sin(pi*x) + c(2)*sin(2*pi*x) + c(3)*sin(3*pi*x);
  assert(max(abs(exact_solution(x, 0, 'coeffs', c) - u0_series)) < 1e-14);
  % each mode decays with its own rate: ratio at t = 0.3
  ratio1 = exact_solution(0.5, 0.3, 'coeffs', [1 0]) ...
         / exact_solution(0.5, 0.0, 'coeffs', [1 0]);
  ratio3 = exact_solution(0.5, 0.3, 'coeffs', [0 0 1]) ...
         / exact_solution(0.5, 0.0, 'coeffs', [0 0 1]);
  assert(abs(ratio1 - exp(-pi^2*0.3)) < 1e-14);
  assert(abs(ratio3 - exp(-9*pi^2*0.3)) < 1e-14);
  % steady linear profile superposition keeps boundary values for all t
  for t = [0 0.2 1.7]
    u = exact_solution(x, t, 'uLeft', 1, 'uRight', 3, 'coeffs', 0.7);
    assert(abs(u(1) - 1) < 1e-12 && abs(u(end) - 3) < 1e-12);
    assert(max(abs(u - (1 + 2*x + 0.7*sin(pi*x)*exp(-pi^2*t)))) < 1e-12);
  end
  fprintf('    series/steady/boundary identities: |err| < 1e-12\n');
end

% ============================================================================
function test_ftcs_vs_analytic()
% (a) FTCS follows the analytic decay at several times within a tolerance
% consistent with its O(dt) + O(dx^2) truncation error.
  nx = 40; tF = 0.2; ns = 720; ce = 180;   % r = alpha*dt/dx^2 = 0.4444
  s = heatwave('nx', nx, 'tFinal', tF, 'nSteps', ns, ...
               'u0', @(x) sin(pi*x), 'captureEvery', ce);
  tol = 6e-4;   % measured max 3.2e-4; O(dt+dx^2)-consistent, see docs
  for tq = [0.05 0.10 0.20]
    k = find(abs(s.t - tq) < 1e-12);
    assert(~isempty(k), 'query time not stored');
    ue = exact_solution(s.x, s.t(k));
    e = max(abs(s.U(k,:) - ue.'));
    fprintf('    t=%.2f  max|FTCS-exact| = %.3e\n', tq, e);
    assert(e < tol);
  end
end

% ============================================================================
function test_implicit_large_dt()
% (b) Backward Euler stays accurate at dt = 0.01, i.e. r = 16 = 32x the
% FTCS stability limit, where explicit FTCS diverges (shown with the
% stability guard explicitly disabled).
  nx = 40; dt = 0.01; tF = 0.2;
  sI = heatwave('nx', nx, 'tFinal', tF, 'dt', dt, 'scheme', 'implicit', ...
                'u0', @(x) sin(pi*x), 'captureEvery', 20);
  assert(abs(sI.info.r - 16) < 1e-9);
  ue = exact_solution(sI.x, tF);
  eI = max(abs(sI.U(end,:) - ue.'));
  fprintf('    implicit at dt=%.3g (r=%.0f): max|err| = %.3e\n', ...
          dt, sI.info.r, eI);
  assert(eI < 2.5e-2);   % measured 1.34e-2; first order in dt: O(dt,dx^2)

  % same dt with FTCS, guard explicitly off, tiny high-mode seed present
  sU = heatwave('nx', nx, 'tFinal', tF, 'dt', dt, 'stabilityGuard', false, ...
                'captureEvery', 20, ...
                'u0', @(x) sin(pi*x) + 0.05*sin(60*pi*x));
  div = max(abs(sU.U(end,:)));
  fprintf('    FTCS at same dt (guard off):  max|u|   = %.3e  (diverged)\n', ...
          div);
  assert(div > 1e6);
end

% ============================================================================
function test_stability_guard()
% (c) The FTCS stability guard rejects dt > dx^2/(2*alpha) with an error;
% it never silently clips.  The exact stability limit itself is accepted.
  nx = 10; tF = 0.06;                 % dx = 0.1, dtMax = dx^2/(2*alpha) = 0.005
  fired = false; id = '';
  try
    heatwave('nx', nx, 'tFinal', tF, 'dt', 0.006, 'u0', @(x) sin(pi*x));
  catch err
    fired = true; id = err.identifier;
  end
  assert(fired, 'guard did not fire for dt > dx^2/(2*alpha)');
  assert(strcmp(id, 'heatwave:unstableStep'));
  s = heatwave('nx', nx, 'tFinal', tF, 'dt', 0.005, 'u0', @(x) sin(pi*x));
  assert(abs(s.info.dt - 0.005) < 1e-15);   % not clipped
  assert(max(abs(s.U(end,:) - exact_solution(s.x, tF).')) < 2e-2);
  fprintf('    guard fired for dt=0.006 > 0.005; dt=0.005 accepted as-is\n');
end

% ============================================================================
function test_convergence_dx()
% (d1) Spatial refinement with r held near 0.45: both FTCS error parts
% scale like dx^2, so the measured order should be ~2.
  nxs = [10 20 40 80 160]; tF = 0.1;
  errs = zeros(size(nxs));
  fprintf('    nx    dx           nSteps  max|err|\n');
  for k = 1:numel(nxs)
    dx = 1/nxs(k);
    ns = ceil(tF / (0.45*dx*dx));
    s = heatwave('nx', nxs(k), 'tFinal', tF, 'nSteps', ns, ...
                 'u0', @(x) sin(pi*x), 'captureEvery', ns);
    assert(s.info.r <= 0.45 + 1e-12);
    errs(k) = max(abs(s.U(end,:) - exact_solution(s.x, tF).'));
    fprintf('    %-5d %-12.3e %-7d %.3e\n', nxs(k), dx, ns, errs(k));
  end
  p = polyfit(log(1./nxs), log(errs), 1);
  p_dx = p(1);
  fprintf('    measured order in dx: %.3f\n', p_dx);
  assert(p_dx > 1.8 && p_dx < 2.2);
end

% ============================================================================
function test_convergence_dt()
% (d2) Temporal refinement at fixed nx.  FTCS has a small, dt-dependent
% error cancellation near r = 1/6, so the clean standard measurement is
% self-convergence: U(dt) against U(dt/2).  Halving dt must halve the
% difference (order ~1).
  nx = 200; dx = 1/nx; tF = 0.05;
  rs = [0.5 0.25 0.125 0.0625 0.03125];
  Ufin = cell(1, numel(rs));
  dtv = zeros(1, numel(rs));
  Eana = zeros(1, numel(rs));
  for k = 1:numel(rs)
    ns = round(tF / (rs(k)*dx*dx));
    s = heatwave('nx', nx, 'tFinal', tF, 'nSteps', ns, ...
                 'u0', @(x) sin(pi*x), 'captureEvery', ns);
    Ufin{k} = s.U(end,:);
    dtv(k) = s.info.dt;
    Eana(k) = max(abs(Ufin{k} - exact_solution(s.x, tF).'));
  end
  npair = numel(rs) - 1;
  pair_errs = zeros(1, npair);
  for k = 1:npair
    pair_errs(k) = max(abs(Ufin{k} - Ufin{k+1}));
  end
  p = polyfit(log(dtv(1:npair)), log(pair_errs), 1);
  p_dt = p(1);
  for k = 1:npair
    fprintf('    dt=%.4e  |U(dt)-U(dt/2)|=%.3e  |U-exact|=%.3e\n', ...
            dtv(k), pair_errs(k), Eana(k));
  end
  fprintf('    measured order in dt (self-convergence): %.3f\n', p_dt);
  assert(p_dt > 0.8 && p_dt < 1.2);
end

% ============================================================================
function test_neumann_mean()
% (e) Insulated ends: both schemes conserve the discrete mean of u to
% roundoff (exact telescoping property of the ghost-cell update), and the
% profile relaxes toward its mean value.
  for sch = {'ftcs', 'implicit'}
    s = heatwave('nx', 60, 'tFinal', 0.5, 'scheme', sch{1}, 'bc', 'neumann', ...
                 'u0', @(x) 1 + 0.5*sin(pi*x), 'nSteps', 4000, ...
                 'captureEvery', 500);
    m0 = mean(s.U(1,:));
    drift = max(abs(mean(s.U, 2) - m0));
    assert(drift < 1e-10);            % measured ~1e-15 (ftcs), ~2e-13 (impl)
    mEnd = mean(s.U(end,:));
    flat = max(abs(s.U(end,:) - mEnd));
    assert(flat < 0.01);              % higher Neumann modes decayed by t=0.5
    fprintf('    %-8s mean drift=%.3e  final spread=%.3e\n', ...
            sch{1}, drift, flat);
  end
end

% ============================================================================
function test_validation()
% (f) Input validation: every malformed input raises a heatwave/exact_
% solution error instead of running.
  bad = { ...
    @() heatwave('nx', 1), ...                                   % nx < 2
    @() heatwave('nx', 10.5), ...                                % non-integer
    @() heatwave('L', 0), ...                                    % L <= 0
    @() heatwave('alpha', -1), ...                               % alpha <= 0
    @() heatwave('tFinal', 0), ...                               % tFinal <= 0
    @() heatwave('nx', 10, 'tFinal', 0.1, 'dt', -0.01), ...      % dt <= 0
    @() heatwave('nx', 10, 'scheme', 'crank'), ...               % bad scheme
    @() heatwave('nx', 10, 'bc', 'robin'), ...                   % bad bc
    @() heatwave('nx', 10, 'banana', 1), ...                     % unknown opt
    @() heatwave('nx', 10, 'dt'), ...                            % odd args
    @() heatwave('nx', 10, 'u0', zeros(1, 10)), ...              % wrong length
    @() heatwave('nx', 10, 'u0', 'text'), ...                    % bad u0 type
    @() heatwave('nx', 10, 'u0', @(x) nan(size(x))), ...         % NaN in u0
    @() heatwave('nx', 10, 'dt', 0.01, 'nSteps', 10), ...        % both steps
    @() heatwave('nx', 10, 'tFinal', 0.1, 'dt', 0.03), ...       % non-dividing
    @() heatwave('nx', 10, 'captureEvery', 0), ...               % bad capture
    @() heatwave('nx', 10, 'stabilityGuard', 'yes'), ...         % not logical
    @() heatwave('nx', 10, 'nx', 20), ...                        % duplicate
    @() exact_solution(0.5, -1), ...                             % t < 0
    @() exact_solution('text', 0.1), ...                         % x not numeric
    @() exact_solution(0.5, 0.1, 'mode', 0), ...                 % mode < 1
    @() exact_solution(0.5, 0.1, 'coeffs', [1 nan]), ...         % NaN coeffs
    @() exact_solution(0.5, 0.1, 'nope', 1)};                    % unknown opt
  for k = 1:numel(bad)
    fired = false; id = '';
    try
      bad{k}();
    catch err
      fired = true; id = err.identifier;
    end
    assert(fired, sprintf('case %d did not raise an error', k));
    assert(strncmp(id, 'heatwave:', 9) || strncmp(id, 'exact_solution:', 15), ...
           sprintf('case %d raised unexpected id "%s"', k, id));
  end
  fprintf('    %d invalid-input cases all rejected\n', numel(bad));
  % valid vector u0 accepted and stored faithfully (dt within the limit)
  s = heatwave('nx', 10, 'tFinal', 0.01, 'nSteps', 5, ...
               'u0', sin(pi*(0:0.1:1)));
  assert(max(abs(s.U(1,:) - sin(pi*(0:0.1:1)))) == 0);
end

function run_results()
%RUN_RESULTS Regenerate every figure and number in docs/results.
%
%   RUN_RESULTS() (from the repository root) reruns the documented cases
%   with the tested runtime (GNU Octave 11.3.0) and writes to
%   docs/results:
%
%     evolution.png        cover: temperature profiles at 13 times
%     evolution_surface.png mesh surface u(x,t) over the stored states
%     convergence_dx.png   max error vs dx (fixed r), measured order ~2
%     convergence_dt.png   self-convergence error vs dt, measured order ~1
%     error_final.png      |numerical - analytic| at the final time
%     frames/frame_*.png   one PNG per computed animation state
%     run_summary.json     machine-readable summary of all measurements
%
%   The animation GIF itself is assembled from the exported frames with
%   tools/make_gif.py (Python + Pillow), because Octave-side animated
%   GIF export is unreliable across builds:
%
%     octave --no-gui --eval "run('run_results.m')"
%     python tools/make_gif.py
%
%   Headless note (Windows): figures are created with
%   'visible'='off' and the gnuplot graphics toolkit, then printed to
%   PNG.  Offscreen Qt PNG export fails on this installation, so the
%   gnuplot path is the documented mechanism.

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
  error('run_results:setup', ...
        'cannot locate src/; run from the repository root.');
end
addpath(srcdir);

outdir = fullfile(fileparts(srcdir), 'docs', 'results');
framedir = fullfile(outdir, 'frames');
if ~exist(outdir, 'dir'); mkdir(outdir); end
if ~exist(framedir, 'dir'); mkdir(framedir); end
% remove stale frames so the GIF never mixes runs
stale = glob(fullfile(framedir, 'frame_*.png'));
for k = 1:numel(stale)
  delete(stale{k});
end

summary = struct();
summary.runtime_name = 'GNU Octave';
summary.runtime_version = version();
summary.generated_localtime = datestr_now();
summary.graphics = 'gnuplot toolkit, visible off, print -dpng';
summary.gif_assembly = ['per-frame PNG export by Octave, assembled with ' ...
                        'tools/make_gif.py (Python Pillow)'];

fprintf('Heatwave results | GNU Octave %s\n', version());
fprintf('output: %s\n', outdir);

% ======================================================================
% Case A: documented reference case, sin(pi x) with zero Dirichlet ends
% ======================================================================
nx = 120; L = 1; alpha = 1; tF = 0.5;
dx = L / nx;
ns = ceil(tF / (0.45 * dx^2 / alpha));       % r slightly below 0.45
ce = floor(ns / 60);                          % ~61 stored states
sol = heatwave('nx', nx, 'tFinal', tF, 'nSteps', ns, 'alpha', alpha, ...
               'u0', @(x) sin(pi*x), 'captureEvery', ce);
fprintf('case A: nx=%d nSteps=%d r=%.4f frames=%d\n', ...
        nx, sol.info.nSteps, sol.info.r, size(sol.U, 1));
summary.case_A = struct('nx', nx, 'nSteps', sol.info.nSteps, ...
  'dt', sol.info.dt, 'r', sol.info.r, 'tFinal', tF, ...
  'frames_stored', size(sol.U, 1));

% ---- cover: profiles at 13 times (blue -> red) ----
sel = round(linspace(1, size(sol.U, 1), 13));
f = new_white_figure();
hold on;
for k = 1:numel(sel)
  c = 0.85 * (k-1) / (numel(sel)-1);
  plot(sol.x, sol.U(sel(k),:), 'color', [c, 0.25*(1-c)+0.1, 0.85*(1-c)+0.1], ...
       'linewidth', 1.6);
end
hold off;
axis([0 L -0.06 1.06]);
xlabel('x / L'); ylabel('u(x,t)');
title(sprintf(['FTCS profiles, nx=%d, r=%.3f (sin \\pi x initial, ' ...
               'Dirichlet ends)'], nx, sol.info.r));
print(f, fullfile(outdir, 'evolution.png'), '-dpng');
close(f);

% ---- mesh surface of the stored evolution ----
sub = 1:3:size(sol.U, 1);
f = new_white_figure();
mesh(sol.x, sol.t(sub), sol.U(sub, :));
xlabel('x / L'); ylabel('t'); zlabel('u(x,t)');
title('u(x,t) from the FTCS run (stored states)');
print(f, fullfile(outdir, 'evolution_surface.png'), '-dpng');
close(f);

% ---- error field at final time ----
uex = exact_solution(sol.x, tF, 'L', L, 'alpha', alpha);
errA = abs(sol.U(end,:) - uex.');
f = new_white_figure();
plot(sol.x, errA, 'b-', 'linewidth', 1.6);
xlabel('x / L'); ylabel('|u_{FTCS} - u_{exact}|');
title(sprintf('FTCS pointwise error at t=%.2f (max %.3e)', ...
              tF, max(errA)));
print(f, fullfile(outdir, 'error_final.png'), '-dpng');
close(f);
summary.case_A.max_error_final_time = max(errA);

% ---- animation frames: every stored state, fixed axes, exact overlay ----
for k = 1:size(sol.U, 1)
  f = new_white_figure();
  plot(sol.x, sol.U(k,:), 'b-;FTCS;', 'linewidth', 2);
  hold on;
  plot(sol.x, exact_solution(sol.x, sol.t(k), 'L', L, 'alpha', alpha), ...
       'r--;analytic;', 'linewidth', 1.4);
  hold off;
  axis([0 L -0.06 1.06]);
  xlabel('x / L'); ylabel('u(x,t)');
  title(sprintf('Heatwave FTCS: nx = %d, r = %.3f, t = %.4f', ...
                nx, sol.info.r, sol.t(k)));
  print(f, fullfile(framedir, sprintf('frame_%03d.png', k)), '-dpng');
  close(f);
end
fprintf('exported %d animation frames\n', size(sol.U, 1));

% ======================================================================
% Case B: convergence studies (same settings as the test suite)
% ======================================================================
% ---- spatial refinement, r held near 0.45 ----
nxs = [10 20 40 80 160]; tFb = 0.1;
errs_dx = zeros(size(nxs)); steps_dx = zeros(size(nxs));
for k = 1:numel(nxs)
  dxb = 1 / nxs(k);
  nsb = ceil(tFb / (0.45 * dxb^2));
  sb = heatwave('nx', nxs(k), 'tFinal', tFb, 'nSteps', nsb, ...
                'u0', @(x) sin(pi*x), 'captureEvery', nsb);
  errs_dx(k) = max(abs(sb.U(end,:) - exact_solution(sb.x, tFb).'));
  steps_dx(k) = nsb;
end
p_dx = polyfit(log(1./nxs), log(errs_dx), 1); p_dx = p_dx(1);
fprintf('convergence dx: order %.3f\n', p_dx);
summary.convergence_dx = struct('nx', nxs, 'dx', 1./nxs, 'nSteps', ...
  steps_dx, 'max_error', errs_dx, 'measured_order', p_dx, ...
  'method', 'fixed r near 0.45, analytic reference, tFinal=0.1');

f = new_white_figure();
loglog(1./nxs, errs_dx, 'bo-;max error;', 'linewidth', 1.6, 'markersize', 6);
hold on;
loglog(1./nxs, errs_dx(1) * ((1./nxs) / (1./nxs(1))).^2, ...
       'k--;slope 2 reference;', 'linewidth', 1.2);
hold off;
xlabel('dx'); ylabel('max |u_{FTCS} - u_{exact}| at t = 0.1');
title(sprintf('spatial convergence (r \\approx 0.45): measured order %.3f', ...
      p_dx));
print(f, fullfile(outdir, 'convergence_dx.png'), '-dpng');
close(f);

% ---- temporal refinement, nx = 200, self-convergence U(dt) vs U(dt/2) ----
nxc = 200; dxc = 1/nxc; tFc = 0.05;
rs = [0.5 0.25 0.125 0.0625 0.03125];
Ufin = cell(1, numel(rs)); dtv = zeros(1, numel(rs));
for k = 1:numel(rs)
  nsc = round(tFc / (rs(k) * dxc^2));
  sc = heatwave('nx', nxc, 'tFinal', tFc, 'nSteps', nsc, ...
                'u0', @(x) sin(pi*x), 'captureEvery', nsc);
  Ufin{k} = sc.U(end,:);
  dtv(k) = sc.info.dt;
end
pair_errs = zeros(1, numel(rs)-1);
for k = 1:numel(pair_errs)
  pair_errs(k) = max(abs(Ufin{k} - Ufin{k+1}));
end
p_dt = polyfit(log(dtv(1:end-1)), log(pair_errs), 1); p_dt = p_dt(1);
fprintf('convergence dt: order %.3f (self-convergence)\n', p_dt);
summary.convergence_dt = struct('nx', nxc, 'dt', dtv(1:end-1), ...
  'self_convergence_error', pair_errs, 'measured_order', p_dt, ...
  'method', 'max|U(dt)-U(dt/2)| at fixed nx=200, tFinal=0.05');

f = new_white_figure();
loglog(dtv(1:end-1), pair_errs, 'bo-;|U(dt) - U(dt/2)|;', ...
       'linewidth', 1.6, 'markersize', 6);
hold on;
loglog(dtv(1:end-1), pair_errs(1) * dtv(1:end-1) / dtv(1), ...
       'k--;slope 1 reference;', 'linewidth', 1.2);
hold off;
xlabel('dt'); ylabel('max |U(dt) - U(dt/2)|');
title(sprintf(['temporal convergence (nx = %d): measured order %.3f ' ...
      '(self-convergence)'], nxc, p_dt));
print(f, fullfile(outdir, 'convergence_dt.png'), '-dpng');
close(f);

% ---- implicit advantage at large dt ----
dt_big = 0.01;                                  % r = 16 = 32x the FTCS limit
si = heatwave('nx', 40, 'tFinal', 0.2, 'dt', dt_big, 'scheme', 'implicit', ...
              'u0', @(x) sin(pi*x), 'captureEvery', 20);
err_impl = max(abs(si.U(end,:) - exact_solution(si.x, 0.2).'));
fprintf('implicit at r=16: max error %.3e\n', err_impl);
summary.implicit_large_dt = struct('nx', 40, 'dt', dt_big, 'r', si.info.r, ...
  'ftcs_limit_dt', si.info.dtMaxStable, 'max_error_vs_analytic', err_impl);

% ---- Neumann mean conservation ----
sn = heatwave('nx', 60, 'tFinal', 0.5, 'bc', 'neumann', ...
              'u0', @(x) 1 + 0.5*sin(pi*x), 'nSteps', 4000, ...
              'captureEvery', 500);
drift = max(abs(mean(sn.U, 2) - mean(sn.U(1,:))));
fprintf('neumann ftcs mean drift: %.3e\n', drift);
summary.neumann_mean_conservation = struct('scheme', 'ftcs', 'nx', 60, ...
  'nSteps', 4000, 'mean_drift', drift);

% ---- write machine-readable summary ----
txt = jsonencode(summary, 'PrettyPrint', true);
fid = fopen(fullfile(outdir, 'run_summary.json'), 'w');
if fid < 0
  error('run_results:io', 'cannot write run_summary.json');
end
fprintf(fid, '%s\n', txt);
fclose(fid);
fprintf('wrote %s\n', fullfile(outdir, 'run_summary.json'));
fprintf('next: python tools/make_gif.py  (assembles animation.gif)\n');
end

% =========================================================================
function f = new_white_figure()
%NEW_WHITE_FIGURE offscreen gnuplot figure, pure white, no grid lines.
  f = figure('visible', 'off');
  graphics_toolkit(f, 'gnuplot');
  set(f, 'color', 'w');
  set(gca, 'color', 'w', 'xgrid', 'off', 'ygrid', 'off', 'zgrid', 'off');
  set(gca, 'fontsize', 11);
end

% =========================================================================
function s = datestr_now()
%DATESTR_NOW local timestamp string, works on Octave and MATLAB.
  s = datestr(now, 'yyyy-mm-dd HH:MM:SS');
end

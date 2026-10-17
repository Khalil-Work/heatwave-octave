function sol = heatwave(varargin)
%SOL = HEATWAVE(Name, Value, ...) Transient 1D heat diffusion solver.
%
%   SOL = HEATWAVE(...) integrates the transient diffusion (heat) equation
%
%       du/dt = ALPHA * d2u/dx2 ,   0 <= x <= L ,  0 <= t <= TFINAL
%
%   on a uniform node-centered grid with two schemes:
%
%     'ftcs'     explicit forward-time centered-space (default).
%                Stability requires dt <= dx^2/(2*ALPHA).  With
%                'stabilityGuard' true (default) an unstable dt raises an
%                error; it is never silently clipped.
%     'implicit' backward Euler with a hand-written Thomas tridiagonal
%                solve.  Unconditionally stable, first order in time.
%
%   Boundary conditions ('bc'):
%
%     'dirichlet'  u(0,t) = ULEFT, u(L,t) = URIGHT (fixed end temperatures)
%     'neumann'    insulated ends, du/dx = 0 at x=0 and x=L (ghost cells).
%                  ULEFT/URIGHT are ignored.  Both schemes conserve the
%                  discrete arithmetic mean of u exactly (to roundoff).
%
%   Name/value options (all optional):
%     'L'             slab length L             (default 1)
%     'nx'            number of grid intervals  (default 50); the grid has
%                     NX+1 nodes, dx = L/NX
%     'alpha'         diffusivity ALPHA > 0     (default 1)
%     'tFinal'        final time TFINAL > 0     (default 0.1)
%     'dt'            time step.  TFINAL/DT must be an integer (within
%                     1e-9 relative); otherwise pass 'nSteps' instead.
%                     If neither 'dt' nor 'nSteps' is given, a safe default
%                     is chosen (FTCS: r = 0.45; implicit: 200 steps).
%     'nSteps'        number of steps; DT = TFINAL/NSTEPS.  Giving both
%                     'dt' and 'nSteps' is an error.
%     'scheme'        'ftcs' (default) or 'implicit'
%     'bc'            'dirichlet' (default) or 'neumann'
%     'u0'            initial condition: function handle of x, or numeric
%                     vector of length NX+1.  Empty uses sin(pi*x/L).
%     'uLeft'         Dirichlet value at x = 0   (default 0)
%     'uRight'        Dirichlet value at x = L   (default 0)
%     'captureEvery'  store the solution every CAPTUREEVERY steps and
%                     always the final step (default 1 = every step)
%     'stabilityGuard' logical (default true).  If true, 'ftcs' errors on
%                     dt > dx^2/(2*ALPHA).  Setting it false is an explicit
%                     opt-out used to demonstrate instability on purpose.
%
%   Returns a struct SOL with fields
%     .x     (NX+1 x 1) node coordinates
%     .t     (nTimes x 1) stored times (includes 0 and TFINAL)
%     .U     (nTimes x NX+1) stored states, row k = u at time SOL.T(k)
%     .info  struct with all solver parameters, the step ratio
%            R = ALPHA*DT/DX^2, and DTMAX = DX^2/(2*ALPHA).
%
%   Core update loops are vectorized.  No toolboxes are used; only core
%   language features, so the code runs on GNU Octave (tested runtime:
%   GNU Octave 11.3.0) and is written to be MATLAB-compatible.
%
%   Example:
%     sol = heatwave('nx',100,'tFinal',0.2,'dt',1e-3,'u0',@(x) sin(pi*x));
%     plot(sol.x, sol.U(end,:));
%
%   See also EXACT_SOLUTION.

% ---------------------------------------------------------------- parsing

opt = struct('L',1, 'nx',50, 'alpha',1, 'tFinal',0.1, 'dt',[], ...
             'nSteps',[], 'scheme','ftcs', 'bc','dirichlet', 'u0',[], ...
             'uLeft',0, 'uRight',0, 'captureEvery',1, ...
             'stabilityGuard',true);

known = {'L','nx','alpha','tFinal','dt','nSteps','scheme','bc','u0', ...
         'uLeft','uRight','captureEvery','stabilityGuard'};
seen  = struct();

nv = numel(varargin);
if mod(nv,2) ~= 0
  error('heatwave:usage', ...
        'heatwave: name/value options must come in pairs.');
end
for k = 1:2:nv
  name = varargin{k};
  if ~ischar(name) || ~any(strcmpi(name, known))
    error('heatwave:invalidOption', ...
          'heatwave: unknown option "%s". Known options: %s.', ...
          tostr(name), strjoin(known', ', '));
  end
  canon = known{find(strcmpi(name, known), 1)};   % canonical spelling
  if isfield(seen, canon)
    error('heatwave:invalidOption', ...
          'heatwave: option ''%s'' given more than once.', canon);
  end
  seen.(canon) = true;
  opt.(canon) = varargin{k+1};
end
opt.scheme = lower(opt.scheme);
opt.bc = lower(opt.bc);

% ------------------------------------------------------------ validation

check_scalar(opt.L,      'L',      true);
check_scalar(opt.alpha,  'alpha',  true);
check_scalar(opt.tFinal, 'tFinal', true);
check_scalar(opt.uLeft,  'uLeft',  false);
check_scalar(opt.uRight, 'uRight', false);

if ~(isnumeric(opt.nx) && isscalar(opt.nx) && isreal(opt.nx) ...
     && opt.nx == round(opt.nx) && opt.nx >= 2)
  error('heatwave:invalidInput', ...
        'heatwave: ''nx'' must be an integer >= 2 (got %s).', ...
        tostr(opt.nx));
end
if ~ischar(opt.scheme) || ~any(strcmp(opt.scheme, {'ftcs','implicit'}))
  error('heatwave:invalidInput', ...
        'heatwave: ''scheme'' must be ''ftcs'' or ''implicit'' (got %s).', ...
        tostr(opt.scheme));
end
if ~ischar(opt.bc) || ~any(strcmp(opt.bc, {'dirichlet','neumann'}))
  error('heatwave:invalidInput', ...
        'heatwave: ''bc'' must be ''dirichlet'' or ''neumann'' (got %s).', ...
        tostr(opt.bc));
end
if ~(isscalar(opt.stabilityGuard) && islogical(opt.stabilityGuard))
  error('heatwave:invalidInput', ...
        'heatwave: ''stabilityGuard'' must be a logical scalar.');
end
if ~(isnumeric(opt.captureEvery) && isscalar(opt.captureEvery) ...
     && isreal(opt.captureEvery) && opt.captureEvery == round(opt.captureEvery) ...
     && opt.captureEvery >= 1)
  error('heatwave:invalidInput', ...
        'heatwave: ''captureEvery'' must be a positive integer.');
end
if ~isempty(opt.dt) && ~isempty(opt.nSteps)
  error('heatwave:invalidInput', ...
        'heatwave: give either ''dt'' or ''nSteps'', not both.');
end

L      = double(opt.L);
nx     = double(opt.nx);
alpha  = double(opt.alpha);
tFinal = double(opt.tFinal);
dx = L / nx;
dtMaxStable = dx^2 / (2*alpha);

% ---- time step selection (never silently modified once chosen) ----
if ~isempty(opt.dt)
  check_scalar(opt.dt, 'dt', true);
  dt = double(opt.dt);
  q = tFinal / dt;
  nSteps = round(q);
  if nSteps < 1 || abs(q - nSteps) > 1e-9 * max(1, abs(q))
    error('heatwave:invalidInput', ...
      ['heatwave: tFinal/dt = %g is not an integer within 1e-9. ' ...
       'Choose dt that divides tFinal exactly, or pass ''nSteps''.'], q);
  end
elseif ~isempty(opt.nSteps)
  if ~(isnumeric(opt.nSteps) && isscalar(opt.nSteps) && isreal(opt.nSteps) ...
       && opt.nSteps == round(opt.nSteps) && opt.nSteps >= 1)
    error('heatwave:invalidInput', ...
          'heatwave: ''nSteps'' must be a positive integer.');
  end
  nSteps = double(opt.nSteps);
  dt = tFinal / nSteps;
else
  % documented defaults: safe step for FTCS (r = 0.45), 200 steps implicit
  if strcmp(opt.scheme, 'ftcs')
    nSteps = max(1, ceil(tFinal / (0.45 * dx^2 / alpha)));
  else
    nSteps = 200;
  end
  dt = tFinal / nSteps;
end
r = alpha * dt / dx^2;

% ---- FTCS stability guard (errors, never clips) ----
guardUsed = false;
if strcmp(opt.scheme, 'ftcs') && opt.stabilityGuard
  guardUsed = true;
  if dt > dtMaxStable * (1 + 1e-9)
    error('heatwave:unstableStep', ...
      ['heatwave: FTCS is unstable for this time step: dt = %.17g exceeds ' ...
       'dx^2/(2*alpha) = %.17g (r = %.17g > 0.5). Reduce dt or nx, or use ' ...
       'scheme=''implicit''.'], dt, dtMaxStable, r);
  end
end

% ---- initial condition ----
x = linspace(0, L, nx+1).';
if isempty(opt.u0)
  u0 = sin(pi * x / L);
elseif isa(opt.u0, 'function_handle')
  u0 = opt.u0(x);
  if ~isnumeric(u0)
    error('heatwave:invalidInput', ...
          'heatwave: ''u0'' function must return numeric values.');
  end
  u0 = double(u0(:));
  if numel(u0) ~= nx+1
    error('heatwave:invalidInput', ...
      'heatwave: ''u0'' function must return one value per grid node.');
  end
elseif isnumeric(opt.u0)
  u0 = double(opt.u0(:));
  if numel(u0) ~= nx+1
    error('heatwave:invalidInput', ...
      'heatwave: ''u0'' vector must have nx+1 = %d entries (got %d).', ...
      nx+1, numel(u0));
  end
else
  error('heatwave:invalidInput', ...
        'heatwave: ''u0'' must be a function handle or a numeric vector.');
end
if ~all(isfinite(u0))
  error('heatwave:invalidInput', ...
        'heatwave: initial condition contains NaN or Inf.');
end
u = u0;

% ---------------------------------------------------------- time stepping

nCap = floor(nSteps / opt.captureEvery) + 1;      % steps 0,ce,2ce,...
if mod(nSteps, opt.captureEvery) ~= 0
  nCap = nCap + 1;                                % plus final partial step
end
tStore = zeros(nCap, 1);
UStore = zeros(nCap, nx+1);
tStore(1) = 0;
UStore(1,:) = u.';
kStore = 1;

switch opt.scheme

  case 'ftcs'
    for n = 1:nSteps
      un = u;
      un(2:nx) = u(2:nx) + r * (u(3:nx+1) - 2*u(2:nx) + u(1:nx-1));
      if strcmp(opt.bc, 'neumann')
        % ghost cells u(0) = u(1), u(nx+2) = u(nx+1): zero-gradient ends
        un(1)    = u(1)    + r * (u(2)    - u(1));
        un(nx+1) = u(nx+1) + r * (u(nx)   - u(nx+1));
      else
        un(1)    = opt.uLeft;
        un(nx+1) = opt.uRight;
      end
      u = un;
      if mod(n, opt.captureEvery) == 0 || n == nSteps
        kStore = kStore + 1;
        tStore(kStore) = n * dt;
        UStore(kStore,:) = u.';
      end
    end

  case 'implicit'  % backward Euler, own Thomas algorithm
    if strcmp(opt.bc, 'dirichlet')
      nUnknown = nx - 1;                    % interior nodes only
      dia = (1 + 2*r) * ones(nUnknown, 1);
      sub = -r * ones(nUnknown, 1);  sub(1) = 0;
      sup = -r * ones(nUnknown, 1);  sup(end) = 0;
    else                                   % neumann: solve all nodes
      nUnknown = nx + 1;
      dia = (1 + 2*r) * ones(nUnknown, 1);
      dia(1) = 1 + r;  dia(end) = 1 + r;    % ghost-cell rows
      sub = -r * ones(nUnknown, 1);  sub(1) = 0;
      sup = -r * ones(nUnknown, 1);  sup(end) = 0;
    end
    for n = 1:nSteps
      if strcmp(opt.bc, 'dirichlet')
        rhs = u(2:nx);
        rhs(1)   = rhs(1)   + r * opt.uLeft;
        rhs(end) = rhs(end) + r * opt.uRight;
        uInt = thomas_solve(sub, dia, sup, rhs);
        u = [opt.uLeft; uInt; opt.uRight];
      else
        u = thomas_solve(sub, dia, sup, u);
      end
      if mod(n, opt.captureEvery) == 0 || n == nSteps
        kStore = kStore + 1;
        tStore(kStore) = n * dt;
        UStore(kStore,:) = u.';
      end
    end

  otherwise
    error('heatwave:invalidInput', 'heatwave: unreachable scheme branch.');
end

% -------------------------------------------------------------- outputs

sol = struct();
sol.x = x;
sol.t = tStore(1:kStore);
sol.U = UStore(1:kStore, :);
sol.info = struct('scheme',opt.scheme, 'bc',opt.bc, 'L',L, 'nx',nx, ...
  'alpha',alpha, 'tFinal',tFinal, 'dt',dt, 'nSteps',nSteps, ...
  'dx',dx, 'r',r, 'dtMaxStable',dtMaxStable, ...
  'captureEvery',opt.captureEvery, 'uLeft',opt.uLeft, ...
  'uRight',opt.uRight, 'stabilityGuard',opt.stabilityGuard, ...
  'guardApplied',guardUsed);
end

% ======================================================================
function w = thomas_solve(sub, dia, sup, rhs)
%THOMAS_SOLVE Solve a tridiagonal system without pivoting (own code).
%   W = THOMAS_SOLVE(SUB, DIA, SUP, RHS) solves
%       SUB(i)*w(i-1) + DIA(i)*w(i) + SUP(i)*w(i+1) = RHS(i)
%   with SUB(1) and SUP(end) ignored.  The matrices built by HEATWAVE are
%   diagonally dominant, so the non-pivoting elimination is stable.
%   Forward sweep eliminates the sub-diagonal; back substitution finishes.
  n = numel(dia);
  b = dia(:);
  c = sup(:);
  d = rhs(:);
  for i = 2:n
    m = sub(i) / b(i-1);
    b(i) = b(i) - m * c(i-1);
    d(i) = d(i) - m * d(i-1);
  end
  w = zeros(n, 1);
  w(n) = d(n) / b(n);
  for i = n-1:-1:1
    w(i) = (d(i) - c(i) * w(i+1)) / b(i);
  end
end

% ======================================================================
function check_scalar(v, nm, must_be_positive)
%CHECK_SCALAR Validate a real finite numeric scalar option.
  ok = isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v);
  if ok && must_be_positive
    ok = v > 0;
  end
  if ~ok
    if must_be_positive
      error('heatwave:invalidInput', ...
            'heatwave: ''%s'' must be a finite positive scalar (got %s).', ...
            nm, tostr(v));
    else
      error('heatwave:invalidInput', ...
            'heatwave: ''%s'' must be a finite scalar (got %s).', ...
            nm, tostr(v));
    end
  end
end

% ======================================================================
function s = tostr(v)
%TOSTR Render an arbitrary value for an error message.
  if ischar(v)
    s = ['''' v ''''];
  elseif isnumeric(v) && isscalar(v)
    s = num2str(v);
  else
    s = class(v);
  end
end

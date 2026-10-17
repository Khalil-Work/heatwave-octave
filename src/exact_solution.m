function u = exact_solution(x, t, varargin)
%U = EXACT_SOLUTION(X, T) Analytic solution of u_t = ALPHA*u_xx on [0,L].
%
%   Reference (documented) case - a single sine mode on a slab with
%   homogeneous Dirichlet ends and initial temperature sin(pi*x/L):
%
%       u(x,t) = sin(pi*x/L) * exp(-ALPHA * pi^2 * t / L^2)
%
%   General decaying-series form.  Separation of variables expands any
%   initial condition that vanishes at x=0 and x=L into sine modes, and
%   each mode decays exponentially with its own rate:
%
%       u(x,t) = S(x) + sum_{k=1}^K  c(k) * sin(k*pi*x/L)
%                                 * exp(-ALPHA * (k*pi/L)^2 * t)
%
%   where S(x) = ULEFT + (URIGHT-ULEFT)*x/L is the steady linear profile
%   associated with non-homogeneous Dirichlet ends (S = 0 when both ends
%   are zero, which the sine series case requires at t = 0).
%
%   Inputs:
%     X   vector of grid points (any orientation; output matches X)
%     T   scalar time (T >= 0)
%
%   Name/value options (all optional):
%     'L'         slab length            (default 1)
%     'alpha'     diffusivity            (default 1)
%     'mode'      sine mode number m     (default 1)
%     'amplitude' amplitude A            (default 1)
%     'coeffs'    coefficient vector c.  If non-empty it overrides
%                 'mode'/'amplitude': c(k) weights mode k = 1..K.
%     'uLeft'     steady-state value at x = 0 (default 0)
%     'uRight'    steady-state value at x = L (default 0)
%
%   Example (documented case, defaults):
%     u = exact_solution(linspace(0,1,101), 0.2);
%
%   See also HEATWAVE.

% ---------------------------------------------------------------- parsing
opt = struct('L',1, 'alpha',1, 'mode',1, 'amplitude',1, 'coeffs',[], ...
             'uLeft',0, 'uRight',0);
known = {'L','alpha','mode','amplitude','coeffs','uLeft','uRight'};

nv = numel(varargin);
if mod(nv,2) ~= 0
  error('exact_solution:usage', ...
        'exact_solution: name/value options must come in pairs.');
end
for k = 1:2:nv
  name = varargin{k};
  if ~ischar(name) || ~any(strcmpi(name, known))
    error('exact_solution:invalidOption', ...
          'exact_solution: unknown option "%s". Known: %s.', ...
          name, strjoin(known', ', '));
  end
  canon = known{find(strcmpi(name, known), 1)};
  if isfield(opt, ['have_' canon]) %#ok<ISFIELD>
    error('exact_solution:invalidOption', ...
          'exact_solution: option ''%s'' given more than once.', canon);
  end
  opt.(['have_' canon]) = true;
  opt.(canon) = varargin{k+1};
end

% ------------------------------------------------------------ validation
if ~(isnumeric(x) && isreal(x) && all(isfinite(x(:))))
  error('exact_solution:invalidInput', ...
        'exact_solution: X must be a finite real numeric vector.');
end
if ~(isnumeric(t) && isscalar(t) && isreal(t) && isfinite(t) && t >= 0)
  error('exact_solution:invalidInput', ...
        'exact_solution: T must be a finite non-negative scalar.');
end
scalar_fields = {'L','alpha','mode','amplitude','uLeft','uRight'};
for k = 1:numel(scalar_fields)
  f = scalar_fields{k};
  v = opt.(f);
  positive = any(strcmp(f, {'L','alpha'})) || strcmp(f,'mode');
  if positive
    if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v > 0 ...
         && ( ~strcmp(f,'mode') || v == round(v) ))
      error('exact_solution:invalidInput', ...
        'exact_solution: ''%s'' must be %s.', f, ...
        'a positive (integer for ''mode'') scalar');
    end
  elseif ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v))
    error('exact_solution:invalidInput', ...
          'exact_solution: ''%s'' must be a finite scalar.', f);
  end
end
if ~isempty(opt.coeffs)
  if ~(isnumeric(opt.coeffs) && isreal(opt.coeffs) ...
       && all(isfinite(opt.coeffs(:))))
    error('exact_solution:invalidInput', ...
          'exact_solution: ''coeffs'' must be a finite real vector.');
  end
end

% -------------------------------------------------------------- evaluate
L = double(opt.L);
alpha = double(opt.alpha);
c = double(opt.coeffs(:));
if isempty(c)
  c = zeros(1, max(1, double(opt.mode)));
  c(double(opt.mode)) = double(opt.amplitude);
end

xv = double(x(:));
K = numel(c);
u = zeros(size(xv));
steady = double(opt.uLeft) + (double(opt.uRight) - double(opt.uLeft)) * xv / L;
for k = 1:K
  if c(k) ~= 0
    u = u + c(k) * sin(k*pi*xv/L) * exp(-alpha * (k*pi/L)^2 * double(t));
  end
end
u = u + steady;

if isvector(x)                    % preserve input orientation
  u = reshape(u, size(x));
end
end

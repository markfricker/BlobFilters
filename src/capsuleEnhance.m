function R = capsuleEnhance(I, params)
% capsuleEnhance  Oriented capsule-kernel ridge/fibre enhancer
%
% USAGE
%   R = capsuleEnhance(I)
%   R = capsuleEnhance(I, params)
%
% DESCRIPTION
%   Enhances elongated bright structures (fibres, neurites, vessels) by
%   convolving the image with a bank of oriented Gaussian capsule kernels
%   across multiple lengths, widths and orientations, then taking the
%   per-pixel maximum response.
%
%   Two kernel modes are available:
%
%   'single'  – each kernel is a plain oriented Gaussian capsule (excitatory
%               centre only).  Zero-mean normalised to suppress DC.
%
%   'doc'     – Difference of Capsules: a narrow excitatory capsule minus a
%               weighted wider inhibitory surround of the same orientation.
%               Provides explicit width selectivity and suppresses slowly
%               varying backgrounds and neighbouring parallel structures,
%               analogous to the role of a Difference of Gaussians (DoG) for
%               blob detection or a Laplacian of Gaussian (LoG) for ridges.
%
% PARAMETERS (fields of the params struct)
%   lengths      – vector of capsule half-lengths in pixels  (default [6 8 12 16])
%   width        – inner (excitatory) capsule width in pixels (default 6)
%   wideWidth    – outer (inhibitory) capsule width for 'doc' mode
%                  (default 2 × width)
%   alpha        – surround weight in 'doc' mode; 0 = no surround,
%                  1 = full subtraction              (default 0.5)
%   orientations – number of evenly spaced angles in [0,180)  (default 8)
%   mode         – 'single' or 'doc'                (default 'single')
%   normalize    – rescale output to [0,1]           (default true)
%
% REFERENCES
%   Lindeberg T. (1998) "Edge detection and ridge detection with automatic
%   scale selection." IJCV 30(2):117-154.
%     → theoretical basis for scale-space ridge operators.
%
%   Frangi A.F. et al. (1998) "Multiscale vessel enhancement filtering."
%   MICCAI, LNCS 1496:130-137.
%     → Hessian-based oriented tubular enhancement; capsule banks pursue
%       a similar goal with explicit kernel shapes rather than Hessian cues.
%
%   Soille P. (2003) "Morphological Image Analysis", 2nd ed., Springer.
%     → oriented structuring elements and granulometry (conceptual basis for
%       length/width scale sweeps).
%
%   Kruger N. & Worgotter F. (2002) "Statistical and deterministic
%   regularities in natural images." PLoS Comput Biol.
%     → oriented elongated filter banks for contour/fibre detection.
%
% NOTES
%   - Input must be a 2-D grayscale image (uint8, uint16, or float).
%   - All computation is performed in single precision.
%   - Requires Image Processing Toolbox for imfilter.
%
% See also: logEnhance, fiberEnhance, fibermetric, imfilter

% -------------------------------------------------------------------------
% defaults
% -------------------------------------------------------------------------
if nargin < 2, params = struct(); end
if ~isfield(params,'lengths'),      params.lengths      = [12 16 20 28 36 40]; end
if ~isfield(params,'width'),        params.width        = 8;           end
if ~isfield(params,'wideWidth'),    params.wideWidth    = params.width * 2.25; end
if ~isfield(params,'alpha'),        params.alpha        = 0.5;         end
if ~isfield(params,'orientations'), params.orientations = 8;           end
if ~isfield(params,'mode'),         params.mode         = 'single';    end
if ~isfield(params,'normalize'),    params.normalize    = true;        end

% -------------------------------------------------------------------------
% input validation
% -------------------------------------------------------------------------
if size(I,3) > 1
    error('capsuleEnhance: expected a 2-D grayscale image.');
end
if ~ismember(lower(params.mode), {'single','doc'})
    error('capsuleEnhance: mode must be ''single'' or ''doc''.');
end
if params.wideWidth <= params.width && strcmpi(params.mode,'doc')
    warning('capsuleEnhance: wideWidth <= width in doc mode — surround may be ineffective.');
end
if any(params.lengths <= 0) || params.width <= 0
    error('capsuleEnhance: lengths and width must be positive.');
end

I = im2single(I);

% -------------------------------------------------------------------------
% pre-build kernel bank (linearised for parfor compatibility)
% -------------------------------------------------------------------------
oris      = linspace(0, 180, params.orientations + 1);
oris(end) = [];
nL        = numel(params.lengths);
nO        = numel(oris);
nKernels  = nL * nO;
kernels   = cell(1, nKernels);

idx = 0;
for ki = 1:nL
    for oi = 1:nO
        idx = idx + 1;
        switch lower(params.mode)
            case 'single'
                kernels{idx} = makeCapsuleKernel( ...
                    params.lengths(ki), params.width, oris(oi));
            case 'doc'
                kernels{idx} = makeDoCapsuleKernel( ...
                    params.lengths(ki), params.width, params.wideWidth, ...
                    oris(oi), params.alpha);
        end
    end
end

% -------------------------------------------------------------------------
% accumulate maximum response across scale/orientation bank
% -------------------------------------------------------------------------
responses = cell(1, nKernels);
if license('test', 'Distrib_Computing_Toolbox')
    parfor ki = 1:nKernels
        responses{ki} = imfilter(I, kernels{ki}, 'replicate', 'conv');
    end
else
    for ki = 1:nKernels
        responses{ki} = imfilter(I, kernels{ki}, 'replicate', 'conv');
    end
end
Rmax = max(cat(3, responses{:}), [], 3);

% -------------------------------------------------------------------------
% normalize
% -------------------------------------------------------------------------
if params.normalize
    Rmax = Rmax - min(Rmax(:));
    rMax = max(Rmax(:));
    if rMax > 0, Rmax = Rmax / rMax; end
end

R = Rmax;
end


% =========================================================================
% kernel helpers  (private, not visible outside this file)
% =========================================================================

function K = makeCapsuleKernel(len, width, angleDeg)
% Single oriented Gaussian capsule kernel (note has to stay double precision for imfilter).
r      = ceil(max(len/2, width/2)) + 6;
[x, y] = meshgrid(-r:r, -r:r);

theta  = -deg2rad(angleDeg);
xr     =  x .* cos(theta) - y .* sin(theta);
yr     =  x .* sin(theta) + y .* cos(theta);

sigmaW = max(width / 2.355, 0.5);
sigmaL = max(len   / 4,     1.0);

K = exp(-(yr .^ 2) ./ (2 * sigmaW ^ 2)) ...
  .* exp(-(xr .^ 2) ./ (2 * sigmaL ^ 2));

K = K - mean(K(:));
K = K / (sum(abs(K(:))) + eps);   % keep double — imfilter requires double kernel
end

% -------------------------------------------------------------------------

function K = makeDoCapsuleKernel(len, wNarrow, wWide, angleDeg, alpha)
% Difference of Capsule kernels — oriented bandpass fibre detector.
%
% The outer (inhibitory) kernel is wider in BOTH axes:
%   - cross-axis : sigWw > sigWn  (suppresses lateral halos)
%   - cap sigma  : sigCapW > sigCapN (suppresses end halos)
%   - body length: outer halfLen is slightly extended so the surround
%                  wraps fully around the flat body and the rounded tips.
% This ensures subtraction occurs uniformly around the entire perimeter
% of the capsule, not just along its sides.

if nargin < 5, alpha = 0.5; end

% Grid sized by the wider kernel
r = ceil(max(len/2, wWide/2)) + 8;
[x, y] = meshgrid(-r:r, -r:r);

theta = -deg2rad(angleDeg);
xr    =  x .* cos(theta) - y .* sin(theta);
yr    =  x .* sin(theta) + y .* cos(theta);

% Sigmas for inner (narrow excitatory) kernel
sigWn   = max(wNarrow / 2.355, 0.5);   % cross-axis
sigCapN = sigWn;                         % cap roundness

% Sigmas for outer (wide inhibitory) kernel
sigWw   = max(wWide / 2.355, 0.5);     % cross-axis — wider laterally
sigCapW = sigWw;                         % cap roundness — wider at ends too

% Inner: flat body + rounded caps using narrow sigmas
halfLenN     = len / 2;
beyondN      = abs(xr) - halfLenN;
K_inner = exp(-(yr.^2) ./ (2*sigWn^2)) ...
        .* exp(-(max(beyondN,0).^2) ./ (2*sigCapN^2));

% Outer: slightly extended body + wider caps — wraps fully around the inner
halfLenW     = len / 2 + (sigWw - sigWn);  % extend body end to match surround
beyondW      = abs(xr) - halfLenW;
K_outer = exp(-(yr.^2) ./ (2*sigWw^2)) ...
        .* exp(-(max(beyondW,0).^2) ./ (2*sigCapW^2));

K = K_inner - alpha .* K_outer;

K = K - mean(K(:));
K = K / (sum(abs(K(:))) + eps);   % keep double — imfilter requires double kernel
end
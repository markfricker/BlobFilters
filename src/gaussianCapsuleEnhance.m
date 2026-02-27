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
if ~isfield(params,'lengths'),      params.lengths      = [6 8 12 16]; end
if ~isfield(params,'width'),        params.width        = 6;           end
if ~isfield(params,'wideWidth'),    params.wideWidth    = params.width * 2; end
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

I = im2single(I);

% -------------------------------------------------------------------------
% pre-build kernel bank
% -------------------------------------------------------------------------
oris     = linspace(0, 180, params.orientations + 1);
oris(end) = [];
nL       = numel(params.lengths);
nO       = numel(oris);
kernels  = cell(nL, nO);

for ki = 1:nL
    for oi = 1:nO
        switch lower(params.mode)
            case 'single'
                kernels{ki,oi} = makeCapsuleKernel( ...
                    params.lengths(ki), params.width, oris(oi));
            case 'doc'
                kernels{ki,oi} = makeDoCapsuleKernel( ...
                    params.lengths(ki), params.width, params.wideWidth, ...
                    oris(oi), single(params.alpha));
        end
    end
end

% -------------------------------------------------------------------------
% accumulate maximum response across scale/orientation bank
% -------------------------------------------------------------------------
Rmax = -inf(size(I), 'single');

for ki = 1:nL
    for oi = 1:nO
        resp = imfilter(I, kernels{ki,oi}, 'replicate', 'conv');
        Rmax = max(Rmax, resp);
    end
end

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
% Single oriented Gaussian capsule kernel (single precision).
r      = ceil(max(len/2, width/2)) + 6;
[x, y] = meshgrid(single(-r:r), single(-r:r));

theta  = single(-deg2rad(angleDeg));
xr     =  x .* cos(theta) - y .* sin(theta);
yr     =  x .* sin(theta) + y .* cos(theta);

sigmaW = single(max(width / 2.355, 0.5));
sigmaL = single(max(len   / 4,     1.0));

K = exp(-(yr .^ 2) ./ (2 * sigmaW ^ 2)) ...
  .* exp(-(xr .^ 2) ./ (2 * sigmaL ^ 2));

K = K - mean(K(:));
K = K / (sum(abs(K(:))) + eps('single'));
end

% -------------------------------------------------------------------------

function K = makeDoCapsuleKernel(len, wNarrow, wWide, angleDeg, alpha)
% Difference of Capsule kernels — oriented bandpass fibre detector.
% Subtracts a weighted wide capsule surround from a narrow excitatory
% centre, both at the same orientation and length.
if nargin < 5, alpha = single(0.5); end

K_inner = makeCapsuleKernel(len, wNarrow, angleDeg);
K_outer = makeCapsuleKernel(len, wWide,   angleDeg);

K = K_inner - alpha .* K_outer;

K = K - mean(K(:));
K = K / (sum(abs(K(:))) + eps('single'));
end
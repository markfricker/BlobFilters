function R = fiberEnhance(I, params)
% fiberEnhance  Multi-scale Hessian-based tubular ridge enhancer
%
% USAGE
%   R = fiberEnhance(I)
%   R = fiberEnhance(I, params)
%
% INPUTS
%   I       - 2-D grayscale image (uint8, uint16, or float); converted to
%             single [0,1] internally.
%   params  - optional struct with fields:
%            .widths    = [6 7 8 9 10]  % cross-sectional widths to probe (px)
%            .multimode = 'builtin'     % aggregation strategy (see below)
%            .normalize = true          % rescale output to [0,1]
%
% OUTPUT
%   R       - enhancement map, single precision, same size as I.
%             Higher values indicate bright, elongated tubular structure.
%
% OVERVIEW
%   Wraps the MATLAB Image Processing Toolbox function fibermetric, which
%   implements a Frangi-style vesselness / tubularity measure based on the
%   eigenvalues of the local image Hessian.
%
%   At each scale w, the Hessian H of the Gaussian-smoothed image is
%   computed (sigma proportional to w).  For a bright tubular ridge, the
%   smaller eigenvalue lambda_1 is large and negative while the larger
%   lambda_2 is near zero; the tubularity score at scale w is:
%
%       S_w = exp(-lambda_1) * (1 - exp(lambda_2 / c))
%
%   where c is an automatic scale factor.  fibermetric aggregates S_w over
%   all requested widths and returns the per-pixel best-scale response.
%
%   MULTIMODE OPTIONS
%   'builtin' (default) — passes all widths to fibermetric at once.
%             MathWorks selects the best-scale response per pixel internally.
%             This is the recommended mode and is fastest.
%   'stack'   — calls fibermetric once per width and takes the pixel-wise
%               maximum.  Mirrors the logEnhance rolling-max pattern and
%               permits inspection of individual scale contributions.
%
%   Width values should span the *short axis* (cross-sectional width) of
%   the target structure, not the long axis.  For 6-8 px wide mitochondria,
%   widths = [6 7 8 9 10] is appropriate.
%
% NOTES
%   - Requires Image Processing Toolbox R2018b or later (fibermetric).
%   - fibermetric internally uses 'bright' object polarity (bright on dark).
%
% REFERENCES
%   Frangi A.F. et al. (1998) Multiscale vessel enhancement filtering.
%   In: Wells W.M. et al. (eds) MICCAI 1998, LNCS 1496, pp. 130-137.
%   Springer, Berlin. https://doi.org/10.1007/BFb0056195
%     -> Original Hessian-based tubularity / vesselness filter.
%
%   MathWorks (2025) fibermetric — Image Processing Toolbox R2025b.
%   https://uk.mathworks.com/help/images/ref/fibermetric.html
%     -> MATLAB implementation used by this wrapper.
%
% EXAMPLE
%   pFib.widths    = [6 7 8 9 10];
%   pFib.multimode = 'builtin';
%   pFib.normalize = true;
%   R = fiberEnhance(I, pFib);
%
% See also: rodGranulometryEnhance, capsuleEnhance, logEnhance

if nargin < 2, params = struct(); end
if ~isfield(params, 'widths'),    params.widths    = [6 7 8 9 10];   end
if ~isfield(params, 'multimode'), params.multimode = 'builtin'; end
if ~isfield(params, 'normalize'), params.normalize = true;        end

% --- input checks -----------------------------------------------------------
if size(I, 3) > 1
    error('fiberEnhance: expected 2-D grayscale image');
end
if ~exist('fibermetric', 'file')
    error('fiberEnhance: fibermetric not found — requires Image Processing Toolbox R2018b+');
end
I = im2single(I);

% --- core enhancement -------------------------------------------------------
switch lower(params.multimode)

    case 'builtin'
        % Pass all widths at once — fibermetric selects the best response
        % per pixel across scales internally. This is the recommended mode.
        R = single(fibermetric(I, params.widths, 'ObjectPolarity', 'bright'));

    case 'stack'
        % Explicit per-width stack, mirroring logEnhance — useful if you
        % want to inspect individual scale responses or weight them.
        R = zeros(size(I), 'single');
        for k = 1:numel(params.widths)
            R = max(R, single(fibermetric(I, params.widths(k), ...
                                   'ObjectPolarity', 'bright')));
        end

    otherwise
        error('fiberEnhance: unknown multimode ''%s''', params.multimode);
end

% --- normalize ---------------------------------------------------------------
if params.normalize
    R = R - min(R(:));
    rMax = max(R(:));
    if rMax > 0, R = R / rMax; end
end
end
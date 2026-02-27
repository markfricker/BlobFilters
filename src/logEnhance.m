function R = logEnhance(I, params)
% logEnhance  Multi-scale Laplacian-of-Gaussian (LoG) enhancement map
%
% R = logEnhance(I, params)
%
% INPUTS
%   I       - grayscale image, expected double in [0,1] (prefiltered externally)
%   params  - optional struct with fields:
%            .sigmas = [2 3 4]        % scales (px) for LoG
%            .normalize = true        % normalize output to [0,1]
%
% OUTPUT
%   R       - enhancement map (higher => more "blob-like" bright centre)
%
% OVERVIEW
%   Computes -LoG at several scales and returns the per-pixel maximum response.
%   Negative sign makes bright blobs have positive response. This is continuous
%   scale-space style enhancement and is robust to noise when combined with
%   mild prefiltering.
%
% REFERENCES
%   - Marr & Hildreth (1980) — Laplacian of Gaussian
%   - Scale-space literature
%
% EXAMPLE
%   R = logEnhance(I, struct('sigmas',[2 3 4 5 6]));
%
% logEnhance  Multi-scale enhancement for elongated bright blobs (mitochondria)

if nargin < 2, params = struct(); end
if ~isfield(params,'sigmas'),    params.sigmas    = [2 3 4 5 6]; end
if ~isfield(params,'normalize'), params.normalize = true;      end

if size(I,3) > 1, error('logEnhance: grayscale input required'); end
I = im2single(I);


rStack = zeros([size(I,1), size(I,2), numel(params.sigmas)],'single');
for k = 1:numel(params.sigmas)
    s    = params.sigmas(k);
    ksz  = 2*ceil(3*s) + 1;
    h    = fspecial('log', ksz, s);
    rStack(:,:,k) = -imfilter(I, h, 'replicate', 'conv');
end
R = max(rStack, [], 3);
R = max(R, 0);  % suppress negative responses (dark blobs)


if params.normalize
    rMin = min(R(:));  R = R - rMin;
    rMax = max(R(:));
    if rMax > 0, R = R / rMax; end
end
end
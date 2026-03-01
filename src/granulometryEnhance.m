function G = granulometryEnhance(I, params)
%
% G = granulometryEnhance(I, params)
%
% INPUTS
%   I       - grayscale image (single or double [0,1]), prefiltered externally
%   params  - optional struct:
%            .sigmas     = [1 2 3 4 8 16 32]   % disk radii in px (width sensitivity)
%            .normalize  = true                 % normalize output to [0,1]
%
% OUTPUT
%   G       - integrated granulometric enhancement map (pattern spectrum sum)
%
% OVERVIEW
%   Computes the granulometric pattern spectrum as the sum of morphological
%   opening residues:
%
%       G = sum_r [ Opening(I, disk(r)) - Opening(I, disk(r+1)) ]
%
%   where r iterates over consecutive pairs in params.sigmas (sorted).
%   Each residue isolates structures at scale r, and summing them produces
%   an enhancement map that highlights objects across all specified scales.
%
%   Note: the final opening (at the largest radius) contributes no residue
%   and is excluded from the sum, consistent with standard granulometry.
%
% NOTES
%   - Sigmas are sorted internally; non-consecutive values are supported but
%     will skip intermediate scales.
%   - For very small objects (< 3 px wide), consider upsampling I before
%     calling this function and mapping results back to original resolution.
%   - The 'orientations' parameter has been removed; disk SEs are isotropic.
%     For anisotropic / line-based granulometry, replace strel('disk',r) with
%     a bank of oriented linear SEs and aggregate over orientations.
%
% REFERENCES
%   Maragos, P. (1989). Pattern spectrum and multiscale shape representation.
%   IEEE Transactions on Pattern Analysis and Machine Intelligence, 11(7),
%   701-716. https://doi.org/10.1109/34.192465
%
%   Serra, J. (1982). Image Analysis and Mathematical Morphology.
%   Academic Press, London.
%

if nargin < 2, params = struct(); end
if ~isfield(params, 'sigmas'),    params.sigmas    = [1 2 3 4 8 16 32]; end
if ~isfield(params, 'normalize'), params.normalize = true;               end

I = im2single(I);
[m, n] = size(I);

% Sort radii so residues are computed between consecutive scales
radii  = sort(params.sigmas(:)', 'ascend');
nRadii = numel(radii);

if any(radii < 1)
    warning('granulometryEnhance:radiusTooSmall', ...
        'Radii < 1 will be clamped to 1; consider larger sigmas or upsampling I first.');
end

if nRadii < 2
    warning('granulometryEnhance:tooFewScales', ...
        'At least 2 sigmas are required to compute residues. Returning zeros.');
    G = zeros(m, n, 'single');
    return
end

G = zeros(m, n, 'single');

% Compute opening at the first (smallest) radius
se_prev  = strel('disk', max(1, round(radii(1))));
open_prev = imopen(I, se_prev);

% Accumulate residues: Opening(r_k) - Opening(r_{k+1})
for k = 2 : nRadii
    se_curr  = strel('disk', max(1, round(radii(k))));
    open_curr = imopen(I, se_curr);

    G = G + (open_prev - open_curr);   % residue at scale r_{k-1}

    % Reuse current opening as previous for next iteration (avoids recompute)
    open_prev = open_curr;
    % se_prev is not reused by imopen, so no need to carry it forward
end

% Normalize to [0, 1]
if params.normalize
    gMin = min(G(:));
    gMax = max(G(:));
    if gMax > gMin
        G = (G - gMin) / (gMax - gMin);
    else
        G = zeros(m, n, 'single');
    end
end

end
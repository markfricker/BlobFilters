function R = rodGranulometryEnhance(I, params)
% rodGranulometryEnhance  Oriented linear granulometry for rod/fibre enhancement
%
% USAGE
%   R = rodGranulometryEnhance(I)
%   R = rodGranulometryEnhance(I, params)
%
% DESCRIPTION
%   Enhances elongated bright structures using a bank of oriented morphological
%   line openings.  At each length L the image is opened with line structuring
%   elements across all specified orientations and the per-pixel maximum is
%   retained:
%
%       O_max(L) = max_θ  imopen(I, strel('line', L, θ))
%
%   The enhancement map is then the sum of pattern-spectrum residues across
%   consecutive length scales (the same formulation used in granulometryEnhance
%   but with line SEs instead of disks):
%
%       R = Σ_k  [ O_max(L_k) − O_max(L_{k+1}) ]
%
%   O_max(L) retains only intensity that belongs to a structure at least L px
%   long along at least one orientation.  The residue isolates structures of
%   length ~L_k.  Unlike disk granulometry, line SEs inherently suppress
%   compact puncta (a round object cannot survive a line opening longer than
%   its diameter), so structureTensorEnhance post-weighting provides less
%   additional benefit here than it does for granulometryEnhance.
%
% PARAMETERS (fields of the params struct)
%   lengths      – vector of line lengths in pixels  (default [8 12 16 20 28 36])
%   orientations – number of evenly spaced angles in [0,180)  (default 8)
%   normalize    – rescale output to [0,1]           (default true)
%
% NOTES
%   - Input must be a 2-D grayscale image (uint8, uint16, or float).
%   - All computation is performed in single precision.
%   - Requires Image Processing Toolbox for imopen and strel.
%   - Parallel execution (parfor) is used automatically if a parallel pool is
%     already open.  Call parpool() before rodGranulometryEnhance to enable
%     this; starting a pool inside the function would cost ~25 s of overhead.
%   - strel('line', L, θ) discretises angle and length to the pixel grid;
%     at orientations near 45 ° the effective SE length may differ by ±1 px.
%     Using ≥ 8 orientations and overlapping lengths averages over this effect.
%   - Lengths are rounded to the nearest integer internally; non-integer inputs
%     are silently coerced.
%
% REFERENCES
%   Maragos P. (1989) Pattern spectrum and multiscale shape representation.
%   IEEE T-PAMI, 11(7):701-716.
%     → theoretical basis for granulometric pattern spectra.
%
%   Soille P. (2003) Morphological Image Analysis, 2nd ed., Springer.
%     → oriented structuring elements and directional granulometry.
%
% See also: granulometryEnhance, capsuleEnhance, structureTensorEnhance

% -------------------------------------------------------------------------
% defaults
% -------------------------------------------------------------------------
if nargin < 2, params = struct(); end
if ~isfield(params, 'lengths'),      params.lengths      = [8 12 16 20 28 36]; end
if ~isfield(params, 'orientations'), params.orientations = 8;                   end
if ~isfield(params, 'normalize'),    params.normalize    = true;                end

% -------------------------------------------------------------------------
% input validation
% -------------------------------------------------------------------------
if size(I, 3) > 1
    error('rodGranulometryEnhance: expected a 2-D grayscale image.');
end
if any(params.lengths < 1)
    error('rodGranulometryEnhance: all lengths must be >= 1.');
end
if params.orientations < 1
    error('rodGranulometryEnhance: orientations must be >= 1.');
end

I = im2single(I);
[m, n] = size(I);

% Sort lengths ascending so residues are computed between consecutive scales
lengths = sort(round(params.lengths(:)'), 'ascend');
nL      = numel(lengths);

if nL < 2
    warning('rodGranulometryEnhance:tooFewLengths', ...
        'At least 2 lengths are required to compute residues. Returning zeros.');
    R = zeros(m, n, 'single');
    return
end

% Orientation bank: [0, 180) exclusive
oris     = linspace(0, 180, params.orientations + 1);
oris(end) = [];
nO       = numel(oris);

% -------------------------------------------------------------------------
% For each length, compute max opening over all orientations.
%
% Use parfor over lengths only if a pool is already open — starting one
% here costs ~25 s and outweighs the benefit for typical bank sizes.
% -------------------------------------------------------------------------
openMax   = cell(1, nL);
useParfor = ~isempty(gcp('nocreate'));

if useParfor
    parfor k = 1:nL
        L  = lengths(k);
        Om = zeros(m, n, 'single');
        for oi = 1:nO
            se = strel('line', L, oris(oi));
            Om = max(Om, imopen(I, se));
        end
        openMax{k} = Om;
    end
else
    for k = 1:nL
        L  = lengths(k);
        Om = zeros(m, n, 'single');
        for oi = 1:nO
            se = strel('line', L, oris(oi));
            Om = max(Om, imopen(I, se));
        end
        openMax{k} = Om;
    end
end

% -------------------------------------------------------------------------
% Accumulate pattern spectrum residues: O_max(L_k) − O_max(L_{k+1})
%
% Each residue highlights structures that survive opening at L_k but not
% at L_{k+1}, i.e., structures of length approximately L_k pixels.
% -------------------------------------------------------------------------
R = zeros(m, n, 'single');
for k = 1 : nL-1
    R = R + (openMax{k} - openMax{k+1});
end

R = max(R, single(0));   % clamp any floating-point negatives

% -------------------------------------------------------------------------
% normalize
% -------------------------------------------------------------------------
if params.normalize
    rMin = min(R(:));
    rMax = max(R(:));
    if rMax > rMin
        R = (R - rMin) / (rMax - rMin);
    else
        R = zeros(m, n, 'single');
    end
end

end

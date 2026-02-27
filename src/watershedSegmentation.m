function labels = watershedSegmentation(I, enhancementMap, seedMask, params)
% watershedSegmentation  Marker-controlled watershed on enhancement map
%
% labels = watershedSegmentation(I, enhancementMap, seedMask, params)
%
% INPUTS
%   I              - original image (for optional intensity-based operations)
%   enhancementMap - continuous enhancement map (0..1)
%   seedMask       - logical mask of seed pixels (markers)
%   params         - struct:
%                    .thrFactor = 0.6     % factor * Otsu for loose mask
%                    .minSize = 20
%
% OUTPUT
%   labels         - labeled segmentation result (integer label image)
%
if nargin < 4, params = struct(); end
if ~isfield(params,'thrFactor'), params.thrFactor = 0.6; end
if ~isfield(params,'minSize'), params.minSize = 12; end

% loose mask from enhancement map
try
    thr = params.thrFactor * graythresh(enhancementMap);  % graythresh gives normalized threshold
catch
    thr = params.thrFactor * 0.5;
end
mask = enhancementMap > thr;
mask = bwareaopen(mask, params.minSize);

% markers - mask seeds to the loose region but keep seeds outside by reinserting
markers = bwlabel(seedMask & mask);
missing = (seedMask & ~mask);
if any(missing(:))
    mn = max(markers(:));
    markers(missing) = (mn+1) * missing(missing);
end

D = bwdist(~mask);
labels = watershed(-D);
% keep only watershed regions where markers exist, then relabel with markers as seeds
labels(~mask) = 0;
% impose markers into watershed result: use marker-based watershed variant
% Simpler: use imimposemin on -D with markers then watershed
minImposed = imimposemin(-D, markers>0);
labels = watershed(minImposed);
labels(~mask) = 0;
% remove small objects and relabel consecutively
labels = bwlabel(labels > 0);
labels = relabelConsecutive(labels);
end

function L2 = relabelConsecutive(L)
    if ~any(L(:)), L2 = L; return; end
    props = regionprops(L, 'PixelIdxList');
    L2 = zeros(size(L));
    for k = 1:numel(props)
        L2(props(k).PixelIdxList) = k;
    end
end
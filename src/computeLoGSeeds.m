function seedMask = computeLoGSeeds(I, params)
% computeLoGSeeds  Compute LoG local maxima seeds (binary mask)
%
% seedMask = computeLoGSeeds(I, params)
%
% INPUTS
%   I      - prefiltered grayscale image (double)
%   params - struct:
%            .sigmas = [2 3 4]
%            .minDistance = 4
%            .thresholdAbs = 0.14
%
% OUTPUT
%   seedMask - logical mask of seed pixels (same size as I)
%
if nargin < 2, params = struct(); end
if ~isfield(params,'sigmas'), params.sigmas = [2 3 4]; end
if ~isfield(params,'minDistance'), params.minDistance = 8; end
if ~isfield(params,'thresholdAbs'), params.thresholdAbs = 0.14; end

R = logEnhance(I, struct('sigmas',params.sigmas,'normalize',true));
peaks = imregionalmax(R);                 % alternative to peak_local_max
% apply threshold and distance pruning
peaks = peaks & (R > params.thresholdAbs);
% enforce min distance by non-max suppression (simple approach)
% Use bwselect on local peaks with morphology (approx)
[rr, cc] = find(peaks);
seedMask = false(size(I));
if ~isempty(rr)
    % greedily add peaks sorted by amplitude
    vals = R(peaks);
    [~, idx] = sort(vals, 'descend');
    used = zeros(size(I));
    for k = 1:numel(idx)
        r = rr(idx(k)); c = cc(idx(k));
        if used(max(1,r-params.minDistance):min(end,r+params.minDistance), max(1,c-params.minDistance):min(end,c+params.minDistance))
            continue;
        else
            seedMask(r,c) = true;
            used(max(1,r-params.minDistance):min(end,r+params.minDistance), max(1,c-params.minDistance):min(end,c+params.minDistance)) = 1;
        end
    end
end
end
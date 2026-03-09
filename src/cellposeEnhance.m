function [BW, L] = cellposeEnhance(I, params)
% cellposeEnhance  Backward-compatible wrapper for cellposeSegment.
%
%   [BW, L] = cellposeEnhance(I)
%   [BW, L] = cellposeEnhance(I, params)
%
%   This function has moved to Segmentation_sandbox/src as cellposeSegment.
%   cellposeEnhance is retained as a thin backward-compatibility shim so
%   that existing code (demos, tests, scripts) continues to work without
%   modification.
%
%   The only difference is output order: cellposeSegment returns [L, BW]
%   (label first), while cellposeEnhance returns [BW, L] (binary first)
%   to preserve the original API.
%
%   New code should call cellposeSegment directly.
%
% See also: cellposeSegment, watershedSegment, refineSegment

if nargin < 2
    [L, BW] = cellposeSegment(I);
else
    [L, BW] = cellposeSegment(I, params);
end
end

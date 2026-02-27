function R = fuseLogPhaseSymEnhance(I, params)
% fuseLogPhaseSymEnhance  Fusion of LoG and PhaseSym enhancement maps
%
% R = fuseLogPhaseSymEnhance(I, params)
%
% INPUTS
%   I      - grayscale image (double)
%   params - struct:
%            .w = 1.0                % weight for PhaseSym map (LoG weight is 1)
%            .logParams = struct()   % passed to logEnhance
%            .psParams = struct()    % passed to phaseSymmetryMinMomentEnhance
%            .normalize = true
%
% OUTPUT
%   R      - fused enhancement map (LoG + w*PhaseSym)
%
% OVERVIEW
%   Computes LoG and PhaseSym maps and returns R = LoG + w*PhaseSym (normalized).
%   In our experiments this fusion with LoG seeds gave improved centroid stability.
%
if nargin < 2, params = struct(); end
if ~isfield(params,'w'), params.w = 1.0; end
if ~isfield(params,'logParams'), params.logParams = struct(); end
if ~isfield(params,'psParams'), params.psParams = struct(); end
if ~isfield(params,'normalize'), params.normalize = true; end

Rlog = logEnhance(I, params.logParams);
[PS, ~] = phaseSymmetryMinMomentEnhance(I, params.psParams);

R = Rlog + params.w * PS;
if params.normalize
    R = R - min(R(:));
    if max(R(:))>0, R = R / max(R(:)); end
end
end
function R = fiberEnhance(I, params)
% fiberEnhance  Multi-scale fibermetric enhancement for mitochondria
%
% R = fiberEnhance(I, params)
%
% INPUTS
%   I       - grayscale image, double or integer class, 2-D
%   params  - optional struct with fields:
%            .widths    = [3 4 5 6]   % object widths in pixels to probe
%            .multimode = 'builtin'   % 'builtin': pass all widths at once to
%                                     %   fibermetric (recommended, MathWorks
%                                     %   implementation selects best scale
%                                     %   per pixel internally)
%                                     % 'stack': compute per-width and take
%                                     %   max (mirrors logEnhance behaviour,
%                                     %   useful for debugging scale responses)
%            .normalize = true        % normalize output to [0,1]
%
% OUTPUT
%   R       - enhancement map in [0,1], higher => more rod-like bright structure
%
% NOTES
%   fibermetric targets tubular/ridge structures via Hessian eigenanalysis
%   (Frangi-style). It is inherently more selective for elongated objects than
%   isotropic LoG, making it better suited to mitochondria.
%
%   Width values should span the *short axis* (width) of the target, not the
%   long axis. For 6-8 px wide mitochondria, widths = [4 5 6 7] is appropriate.
%   Erring slightly narrow (3-6) avoids responding to larger organelles.
%
%   Requires Image Processing Toolbox R2018b or later.
%
% EXAMPLE
%   R = fiberEnhance(I, struct('widths', [3 4 5 6], 'multimode', 'builtin'));

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
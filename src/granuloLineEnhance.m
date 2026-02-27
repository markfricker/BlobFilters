function G = granuloLineEnhance(I, params)
% granuloLineEnhance  Disk + Line granulometry enhancement (granulo-bowler hybrid)
%
% G = granuloLineEnhance(I, params)
%
% INPUTS
%   I       - grayscale image (double [0,1]) prefiltered externally
%   params  - optional struct:
%            .diskRadii = [1 2 3]            % disk radii in px (width sensitivity)
%            .lineLengths = [12 16]          % line lengths in px (elongation)
%            .orientations = 8               % number of orientations
%            .useReconstruction = false     % use opening-by-reconstruction
%            .normalize = true
%
% OUTPUT
%   G       - combined granulometric enhancement map
%
% OVERVIEW
%   Computes disk residues: Opening(disk(r)) - Opening(disk(r+1)) for each radius r,
%   sums them -> Sd. Computes line residues for several lengths and orientations,
%   takes orientation max per length -> Sl. Then G = Sd .* Sl (pixelwise product).
%
% NOTES
%   - This is morphological and discrete; performance depends on structuring elements.
%   - For small objects (few pixels wide) consider upsampling before calling this
%     function and mapping seeds/labels back to original resolution.
%
% REFERENCES
%   - Fricker & Obara (granulometry)
%   - Bowler-hat morphological filters
%
if nargin < 2, params = struct(); end
if ~isfield(params,'diskRadii'), params.diskRadii = [1 2 3 4 8 16 32]; end
if ~isfield(params,'lineLengths'), params.lineLengths = [8 12 16]; end
if ~isfield(params,'orientations'), params.orientations = 16; end
if ~isfield(params,'useReconstruction'), params.useReconstruction = false; end
if ~isfield(params,'normalize'), params.normalize = true; end

I = double(I);
[m,n] = size(I);

% Disk residues
Sd = zeros(m,n);
for r = params.diskRadii
    se1 = strel('disk', max(1,round(r)));
    % se2 = strel('disk', max(1,round(r))+1);
    if params.useReconstruction
        marker1 = imerode(I, se1);
        open1 = imreconstruct(marker1, I);
        % marker2 = imerode(I, se2);
        % open2 = imreconstruct(marker2, I);
    else
        open1 = imopen(I, se1);
        % open2 = imopen(I, se2);
    end
    Sd = Sd + open1;%(open1 - open2);
end

% % Line residues (multi-orientation)
% Sl = zeros(m,n);
% oris = linspace(0,180,params.orientations+1);
% oris(end) = [];
% for L = params.lineLengths
%     orientResponses = zeros(m,n, numel(oris));
%     for k = 1:numel(oris)
%         ang = oris(k);
%         % Build a thin line strel by creating a rotated line in a small image
%         % Use imrotate on a straight line mask for subpixel-ish SE
%         len = round(L);
%         lineMask = zeros(len, len);
%         cx = ceil(len/2); cy = cx;
%         for x = 1:len, lineMask(cx,x) = 1; end
%         lineMask = imrotate(lineMask, ang, 'bicubic','crop') > 0.5;
%         se1 = strel(lineMask);
%         % slightly longer for residue
%         len2 = max(1, round(L+2));
%         lineMask2 = zeros(len2, len2);
%         cx = ceil(len2/2);
%         for x = 1:len2, lineMask2(cx,x) = 1; end
%         lineMask2 = imrotate(lineMask2, ang, 'bicubic','crop') > 0.5;
%         se2 = strel(lineMask2);
%         if params.useReconstruction
%             m1 = imerode(I, se1); o1 = imreconstruct(m1, I);
%             m2 = imerode(I, se2); o2 = imreconstruct(m2, I);
%         else
%             o1 = imopen(I, se1);
%             o2 = imopen(I, se2);
%         end
%         orientResponses(:,:,k) = (o1 - o2);
%     end
%     % max across orientations to emphasize elongation any angle
%     Sl = Sl + max(orientResponses, [], 3);
% end

G = Sd;% .* Sl;
if params.normalize
    G = G - min(G(:));
    if max(G(:))>0, G = G / max(G(:)); end
end
end
function [PSmap, MMmap] = phaseSymmetryMinMomentEnhance(I, params)
% phaseSymmetryMinMomentEnhance
% Approximate Phase Symmetry and Minimum Moment maps using MATLAB gabor()
%
% NOTE:
%   MATLAB gabor uses WAVELENGTH (>=2 pixels), not frequency.
%   This version converts frequency to wavelength properly.
%
% REFERENCES:
%   P. Kovesi – Phase Congruency and Minimum Moment
%

if nargin < 2, params = struct(); end
if ~isfield(params,'wavelengths'), params.wavelengths = [8,16,32,64]; end
if ~isfield(params,'orientations'), params.orientations = 12; end
if ~isfield(params,'tophatRadius'), params.tophatRadius = 16; end
if ~isfield(params,'normalize'), params.normalize = true; end

I = double(I);
oris = linspace(0, pi, params.orientations+1); oris(end) = [];
nOr = numel(oris);

sumEven = zeros(size(I));
sumOdd  = zeros(size(I));
sumAmpPerOr = zeros([size(I), nOr]);

for w = params.wavelengths
    lambda = max(2, round(w));
    for k = 1:nOr
        theta = rad2deg(oris(k));
        g = gabor(lambda, theta);
        R = imgaborfilt(I, g);
        even = real(R);
        odd  = imag(R);
        amp  = sqrt(even.^2 + odd.^2);
        sumEven = sumEven + even;
        sumOdd  = sumOdd  + odd;
        sumAmpPerOr(:,:,k) = sumAmpPerOr(:,:,k) + amp;
    end
end

sumAmp = sum(sumAmpPerOr,3) + eps;
PSmap = (abs(sumEven) - abs(sumOdd)) ./ sumAmp;
PSmap = PSmap - min(PSmap(:));
if params.normalize && max(PSmap(:))>0
    PSmap = PSmap / max(PSmap(:));
end

% Minimum moment tensor
cosT = cos(oris); sinT = sin(oris);
M11 = zeros(size(I)); M22 = zeros(size(I)); M12 = zeros(size(I));
for k = 1:nOr
    E = sumAmpPerOr(:,:,k);
    M11 = M11 + E .* (cosT(k)^2);
    M22 = M22 + E .* (sinT(k)^2);
    M12 = M12 + E .* (cosT(k)*sinT(k));
end

traceM = M11 + M22;
detM = M11 .* M22 - M12 .* M12;
disc = traceM.^2 - 4*detM;
disc(disc < 0) = 0;
eig2 = 0.5 .* (traceM - sqrt(disc));

MMmap = eig2 - min(eig2(:));
if params.normalize && max(MMmap(:))>0
    MMmap = MMmap / max(MMmap(:));
end

% Background suppression
if params.tophatRadius > 0
    se = strel('disk', params.tophatRadius);
    wh = imtophat(I, se);
    if max(wh(:))>0
        wh = (wh - min(wh(:))) / max(wh(:));
        PSmap = PSmap .* wh;
        MMmap = MMmap .* wh;
    end
end

if params.normalize
    PSmap = PSmap - min(PSmap(:));
    if max(PSmap(:))>0, PSmap = PSmap / max(PSmap(:)); end
    MMmap = MMmap - min(MMmap(:));
    if max(MMmap(:))>0, MMmap = MMmap / max(MMmap(:)); end
end
end
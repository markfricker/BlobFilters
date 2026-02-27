function R = hessianBlobnessEnhance(I, params)
% hessianBlobnessEnhance  Hessian-based multi-scale blobness enhancement
%
% R = hessianBlobnessEnhance(I, params)
%
% INPUTS
%   I      - grayscale image (double)
%   params - struct:
%            .scales = [1 1.5 2.0]   % Gaussian scales for second derivatives
%            .method = 'det'         % 'det' or 'minE' (min eigenvalue)
%            .normalize = true
%
% OUTPUT
%   R      - normalized enhancement map
%
% OVERVIEW
%   Computes the Hessian matrix (Ixx, Ixy, Iyy) at multiple scales (via Gaussian
%   smoothing + derivative kernels), computes either determinant or min eigenvalue,
%   and returns the pixelwise maximum across scales.
%
% REFERENCE
%   - Frangi, Sato (vesselness) and multi-scale Hessian analysis
if nargin < 2, params = struct(); end
if ~isfield(params,'scales'), params.scales = [6 8 12]; end
if ~isfield(params,'method'), params.method = 'det'; end
if ~isfield(params,'normalize'), params.normalize = true; end

I = double(I);
Rstack = -inf([size(I), numel(params.scales)]);

for k = 1:numel(params.scales)
    s = params.scales(k);
    % compute second derivatives via Gaussian derivatives
    % Create derivative filters
    hxx = fspecial('gaussian', max(3,ceil(6*s)), s);
    % Derivatives via imfilter using finite differences on smoothed image
    Is = imgaussfilt(I, s);
    [Ix, Iy] = gradient(Is);
    [Ixx, Ixy] = gradient(Ix);
    [~, Iyy] = gradient(Iy);
    % Alternative more accurate second derivatives could be used
    % Compute eigenvalues of Hessian matrix per pixel
    tmp1 = (Ixx + Iyy) / 2;
    tmp2 = sqrt( ((Ixx - Iyy)/2).^2 + Ixy.^2 );
    lambda1 = tmp1 + tmp2;
    lambda2 = tmp1 - tmp2;
    switch lower(params.method)
        case 'det'
            val = lambda1 .* lambda2;
        case 'mine'
            val = lambda2; % assuming |lambda1|>=|lambda2|
        otherwise
            val = lambda1 .* lambda2;
    end
    Rstack(:,:,k) = val;
end

R = max(Rstack, [], 3);
if params.normalize
    R = R - min(R(:));
    if max(R(:))>0, R = R / max(R(:)); end
end
end
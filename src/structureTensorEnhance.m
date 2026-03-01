function [C, theta] = structureTensorEnhance(I, params)
% structureTensorEnhance  Local structure tensor coherence map
%
% USAGE
%   C          = structureTensorEnhance(I)
%   C          = structureTensorEnhance(I, params)
%   [C, theta] = structureTensorEnhance(I, params)
%
% DESCRIPTION
%   Computes the coherence index C of the local structure tensor at each
%   pixel.  C lies in [0,1] and characterises local image geometry:
%
%     C → 1   elongated linear structure (rods, fibres, network segments)
%     C → 0   isotropic region (noise, puncta, junctions, flat background)
%
%   The coherence index is derived from the eigenvalues λ₁ ≥ λ₂ of the
%   smoothed gradient outer-product tensor J:
%
%       C = ( (λ₁ − λ₂) / (λ₁ + λ₂) )²
%
%   C is intended to be used as a multiplicative weight applied to an
%   existing enhancement map AFTER it has been computed:
%
%       R_rod = capsuleEnhance(I, pCap) .* structureTensorEnhance(I, pST);
%
%   This suppresses responses in noise and puncta regions while preserving
%   rod and network responses.  It does not modify the input image passed
%   to the enhancers.
%
% PARAMETERS (fields of the params struct)
%   sigmaGrad  – gradient scale (px): Gaussian smoothing applied before
%                differentiation, suppresses pixel noise.
%                Typical value: 1–2 px.        (default 1.5)
%   sigmaInt   – integration scale (px): Gaussian window used to smooth
%                the tensor components, sets the neighbourhood over which
%                orientation is averaged.  Should match the structure
%                cross-section width (~half-width of target).
%                For 8 px wide mitochondria: 4–6 px.  (default 5)
%   normalize  – rescale C to [0,1]            (default true)
%
% OUTPUTS
%   C     – coherence map, single precision, same size as I
%   theta – dominant orientation map (radians, along-fibre direction),
%           in [-π/2, π/2].  Only computed if requested as second output.
%
% NOTES
%   - Input must be a 2-D grayscale image (uint8, uint16, or float).
%   - All computation is performed in single precision.
%   - Requires Image Processing Toolbox for imfilter.
%   - The coherence index is dimensionless and illumination-invariant:
%     scaling I by a constant does not change C.
%   - For a purely isotropic blob (round punctum) surrounded by dark
%     background, gradients point radially outward → λ₁ ≈ λ₂ → C ≈ 0.
%   - For a long rod, gradients point consistently perpendicular to the
%     long axis → λ₁ >> λ₂ → C → 1.
%
% REFERENCES
%   Weickert J. (1998) Anisotropic Diffusion in Image Processing. Teubner.
%     → original structure tensor formulation and coherence index.
%
%   Bigun J. & Granlund G.H. (1987) Optimal orientation detection of
%   linear symmetry. ICCV, pp. 433-438.
%     → early use of the structure tensor for orientation estimation.
%
% See also: capsuleEnhance, logEnhance, fiberEnhance, granulometryEnhance

% -------------------------------------------------------------------------
% defaults
% -------------------------------------------------------------------------
if nargin < 2, params = struct(); end
if ~isfield(params, 'sigmaGrad'), params.sigmaGrad = 1.5; end
if ~isfield(params, 'sigmaInt'),  params.sigmaInt  = 5;   end
if ~isfield(params, 'normalize'), params.normalize = true; end

% -------------------------------------------------------------------------
% input validation
% -------------------------------------------------------------------------
if size(I, 3) > 1
    error('structureTensorEnhance: expected a 2-D grayscale image.');
end
if params.sigmaGrad <= 0 || params.sigmaInt <= 0
    error('structureTensorEnhance: sigmaGrad and sigmaInt must be positive.');
end

I = im2single(I);

% -------------------------------------------------------------------------
% step 1: smooth at gradient scale then differentiate
%
% Smooth with a separable Gaussian at scale sigmaGrad, then compute
% central-difference gradients.  This is equivalent to convolving with
% the derivative of a Gaussian but avoids implementing a separate
% derivative kernel.
% -------------------------------------------------------------------------
sg   = params.sigmaGrad;
kszG = 2*ceil(3*sg) + 1;
gG   = fspecial('gaussian', [1, kszG], sg);          % double kernel
Ig   = imfilter(imfilter(I, gG, 'replicate'), gG', 'replicate');

[Ix, Iy] = gradient(Ig);   % central differences; output is single

% -------------------------------------------------------------------------
% step 2: form and smooth structure tensor components at integration scale
%
% J = Kρ * (∇I ∇Iᵀ)  where Kρ denotes Gaussian smoothing at scale sigmaInt
%
% Components:
%   J11 = Kρ * Ix²
%   J12 = Kρ * (Ix·Iy)
%   J22 = Kρ * Iy²
% -------------------------------------------------------------------------
si   = params.sigmaInt;
kszI = 2*ceil(3*si) + 1;
gI   = fspecial('gaussian', [1, kszI], si);          % double kernel

J11 = imfilter(imfilter(Ix.*Ix, gI, 'replicate'), gI', 'replicate');
J12 = imfilter(imfilter(Ix.*Iy, gI, 'replicate'), gI', 'replicate');
J22 = imfilter(imfilter(Iy.*Iy, gI, 'replicate'), gI', 'replicate');

% -------------------------------------------------------------------------
% step 3: compute eigenvalues analytically
%
% For symmetric 2×2 matrix [[J11, J12]; [J12, J22]]:
%   discriminant = sqrt((J11-J22)² + 4·J12²)   ← numerically stable form
%   λ₁ = (trace + discriminant) / 2             ← larger eigenvalue
%   λ₂ = (trace - discriminant) / 2             ← smaller eigenvalue
% -------------------------------------------------------------------------
disc = sqrt((J11 - J22).^2 + 4*J12.^2);   % always ≥ 0
lam1 = (J11 + J22 + disc) / 2;            % larger  (across fibre)
lam2 = (J11 + J22 - disc) / 2;            % smaller (along  fibre)

% -------------------------------------------------------------------------
% step 4: coherence index
%
%   C = ( (λ₁ − λ₂) / (λ₁ + λ₂) )²
%
% C is in [0,1].  The small epsilon in the denominator guards against
% dividing by zero in flat background regions where λ₁ + λ₂ ≈ 0.
% -------------------------------------------------------------------------
denom = lam1 + lam2 + single(1e-10);
C     = ((lam1 - lam2) ./ denom) .^ 2;
C     = max(C, single(0));   % clamp any floating-point negatives

% -------------------------------------------------------------------------
% step 5: orientation map (optional second output)
%
% The along-fibre direction is the eigenvector of the SMALLER eigenvalue
% (λ₂).  For a 2×2 symmetric matrix the angle of the eigenvector of λ₁
% (across-fibre) is:
%
%   θ_across = 0.5 · atan2(2·J12, J11 − J22)
%
% The along-fibre direction is perpendicular:
%
%   θ_fibre = θ_across + π/2,  wrapped to [−π/2, π/2]
% -------------------------------------------------------------------------
if nargout > 1
    theta_across = single(0.5) .* atan2(2*J12, J11 - J22);
    theta_fibre  = theta_across + single(pi/2);
    % Wrap to [-π/2, π/2]:  mod(θ + π/2, π) − π/2
    theta = mod(theta_fibre + single(pi/2), single(pi)) - single(pi/2);
end

% -------------------------------------------------------------------------
% normalize
% -------------------------------------------------------------------------
if params.normalize
    cMax = max(C(:));
    if cMax > 0, C = C / cMax; end
end

end

function Ism = orientedGaussSmooth(I, params)
% orientedGaussSmooth  Orientation-adaptive Gaussian smoothing
%
% USAGE
%   Ism = orientedGaussSmooth(I)
%   Ism = orientedGaussSmooth(I, params)
%
% DESCRIPTION
%   Smooths the image with a Gaussian that is elongated along the local
%   fibre/rod orientation estimated from the structure tensor.  In flat or
%   isotropic regions (noise, puncta) the filter falls back to isotropic
%   Gaussian smoothing; in linear regions it extends strongly along the
%   fibre axis and barely across it.
%
%   This is a single-pass, non-iterative approximation to Coherence-Enhancing
%   Diffusion (CED).  Multiple passes (feed output back as input) produce
%   stronger gap-filling at the cost of additional compute time.
%
%   ALGORITHM
%   1. Estimate local orientation θ and coherence C from the structure tensor.
%   2. Pre-convolve I with N oriented Gaussians at discrete angles θ_k.
%   3. Blend the convolved images per-pixel using cosine-squared weights
%      scaled by the coherence index:
%
%        w_k(x,y) = C(x,y)·cos²(θ(x,y) − θ_k) + (1−C(x,y))/N
%        Ism      = Σ_k w_k · G_k(I) / Σ_k w_k
%
%   Because Σ_k cos²(θ − θ_k) = N/2 exactly for N uniformly-spaced θ_k,
%   the denominator simplifies to 1 + C·(N/2 − 1) (no per-pixel loop).
%   The (1−C)/N isotropic component ensures the output is well-defined and
%   noise-smoothed in flat background regions.
%
% PARAMETERS (fields of the params struct)
%   sigmaAlong   – along-fibre smoothing scale (px)   (default 4)
%   sigmaAcross  – across-fibre smoothing scale (px)  (default 1.5)
%   orientations – discrete orientations in [0,180)   (default 8)
%   sigmaGrad    – structure tensor gradient scale (px)(default 1.5)
%   sigmaInt     – structure tensor integration scale  (default 5)
%
% OUTPUTS
%   Ism  – smoothed image, single precision, same size as I.
%          Pixel values are a weighted blend of oriented smoothings so
%          intensities are preserved (no amplification).
%
% NOTES
%   - Input must be a 2-D grayscale image (uint8, uint16, or float).
%   - All computation is performed in single precision.
%   - Requires Image Processing Toolbox for imfilter.
%   - imfilter requires double kernels; kernels are built in double and
%     passed to imfilter, which returns a single output matching I.
%   - Parallel execution (parfor) is used automatically if a parallel pool
%     is already open.
%
% REFERENCES
%   Weickert J. (1999) Coherence-enhancing diffusion filtering.
%   IJCV 31(2-3):111-127.
%     → original CED; this function approximates one PDE iteration.
%
%   Freeman W.T. & Adelson E.H. (1991) The design and use of steerable
%   filters. IEEE T-PAMI 13(9):891-906.
%     → theoretical basis for blending discrete orientation responses.
%
% See also: structureTensorEnhance, imdiffusefilt, imguidedfilter

% -------------------------------------------------------------------------
% defaults
% -------------------------------------------------------------------------
if nargin < 2, params = struct(); end
if ~isfield(params, 'sigmaAlong'),   params.sigmaAlong   = 4;   end
if ~isfield(params, 'sigmaAcross'),  params.sigmaAcross  = 1.5; end
if ~isfield(params, 'orientations'), params.orientations = 8;   end
if ~isfield(params, 'sigmaGrad'),    params.sigmaGrad    = 1.5; end
if ~isfield(params, 'sigmaInt'),     params.sigmaInt     = 5;   end

% -------------------------------------------------------------------------
% input validation
% -------------------------------------------------------------------------
if size(I, 3) > 1
    error('orientedGaussSmooth: expected a 2-D grayscale image.');
end
if params.sigmaAlong <= 0 || params.sigmaAcross <= 0
    error('orientedGaussSmooth: sigmaAlong and sigmaAcross must be positive.');
end
if params.sigmaAlong < params.sigmaAcross
    warning('orientedGaussSmooth:sigmaOrder', ...
        'sigmaAlong < sigmaAcross — the filter will not be elongated along fibres.');
end

I = im2single(I);
[m, n] = size(I);

% -------------------------------------------------------------------------
% step 1: structure tensor → local coherence C and along-fibre angle theta
%
% Identical derivation to structureTensorEnhance; inline to avoid a
% function dependency and to share gradient computations.
% -------------------------------------------------------------------------
sg   = params.sigmaGrad;
kszG = 2*ceil(3*sg) + 1;
gG   = fspecial('gaussian', [1, kszG], sg);          % 1-D, double
Ig   = imfilter(imfilter(I, gG, 'replicate'), gG', 'replicate');

[Ix, Iy] = gradient(Ig);   % central differences, single output

si   = params.sigmaInt;
kszI = 2*ceil(3*si) + 1;
gI   = fspecial('gaussian', [1, kszI], si);           % 1-D, double

J11 = imfilter(imfilter(Ix.*Ix, gI, 'replicate'), gI', 'replicate');
J12 = imfilter(imfilter(Ix.*Iy, gI, 'replicate'), gI', 'replicate');
J22 = imfilter(imfilter(Iy.*Iy, gI, 'replicate'), gI', 'replicate');

disc = sqrt((J11 - J22).^2 + 4*J12.^2);
lam1 = (J11 + J22 + disc) / 2;
lam2 = (J11 + J22 - disc) / 2;

% Coherence index C ∈ [0,1]
C = ((lam1 - lam2) ./ (lam1 + lam2 + single(1e-10))) .^ 2;
C = max(C, single(0));

% Along-fibre orientation theta ∈ [-π/2, π/2]
theta_across = single(0.5) .* atan2(2*J12, J11 - J22);
theta_fibre  = theta_across + single(pi/2);
theta = mod(theta_fibre + single(pi/2), single(pi)) - single(pi/2);

% -------------------------------------------------------------------------
% step 2: build bank of oriented Gaussian kernels and pre-convolve
%
% Each kernel is elongated along its orientation axis:
%   - long axis  (xr): sigma = sigmaAlong  (along the fibre)
%   - short axis (yr): sigma = sigmaAcross (across the fibre)
%
% Kernel at 0° is elongated horizontally; at 90°, vertically.
% Kernels are kept double — imfilter requires a double filter argument.
% -------------------------------------------------------------------------
N        = params.orientations;
oris_deg = linspace(0, 180, N + 1);
oris_deg(end) = [];                       % [0, 180) exclusive, length N
oris_rad = single(oris_deg * pi / 180);   % radians, for blending step

sA = params.sigmaAlong;
sX = params.sigmaAcross;
r  = ceil(3 * sA) + 1;
[xg, yg] = meshgrid(-r:r, -r:r);         % double grid

Gk_cell   = cell(1, N);
useParfor = ~isempty(gcp('nocreate'));

if useParfor
    parfor k = 1:N
        th = -oris_deg(k) * pi / 180;           % rotation (negative = CCW in image coords)
        xr =  xg * cos(th) - yg * sin(th);      % along-kernel axis
        yr =  xg * sin(th) + yg * cos(th);      % across-kernel axis
        K  = exp(-xr.^2 / (2*sA^2) - yr.^2 / (2*sX^2));
        K  = K / sum(K(:));                      % unit-sum: smoothing, not detection
        Gk_cell{k} = imfilter(I, K, 'replicate');
    end
else
    for k = 1:N
        th = -oris_deg(k) * pi / 180;
        xr =  xg * cos(th) - yg * sin(th);
        yr =  xg * sin(th) + yg * cos(th);
        K  = exp(-xr.^2 / (2*sA^2) - yr.^2 / (2*sX^2));
        K  = K / sum(K(:));
        Gk_cell{k} = imfilter(I, K, 'replicate');
    end
end

% -------------------------------------------------------------------------
% step 3: coherence-weighted cosine² blending
%
%   w_k  = C .* cos²(theta − theta_k) + (1−C) / N
%
%   Σ_k cos²(theta − theta_k) = N/2  (exactly, for uniform spacing)
%   ⟹  Σ_k w_k = 1 + C .* (N/2 − 1)   (pixel-wise scalar, no loop needed)
%
%   High C → weight concentrated near local orientation (steered smoothing)
%   Low  C → equal weights for all orientations (isotropic smoothing)
% -------------------------------------------------------------------------
Inumer = zeros(m, n, 'single');
for k = 1:N
    wk     = C .* (cos(theta - oris_rad(k)) .^ 2) + (single(1) - C) ./ single(N);
    Inumer = Inumer + wk .* Gk_cell{k};
end
Idenom = single(1) + C .* (single(N)/2 - single(1));

Ism = Inumer ./ Idenom;

end

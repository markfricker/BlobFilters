function I = makeMitoTestImage(sz, seed)
% makeMitoTestImage  Synthetic confocal-like mitochondria test image
%
% USAGE
%   I = makeMitoTestImage()           % 512x512, fixed seed
%   I = makeMitoTestImage(sz, seed)
%
% DESCRIPTION
%   Generates a single-precision grayscale image containing structures
%   matched to real confocal mitochondria dimensions:
%
%   Isolated structures (sparse region, top-left quadrant)
%     - Puncta        : round mitochondria, diameters 6-10 px
%     - Short capsules: early-elongation forms, lengths 12-20 px, width ~8 px
%     - Long capsules : elongated rods, lengths 20-40 px, width ~8 px
%
%   Clustered structures (centre region)
%     - Dense packing of mixed short capsules and puncta with slight
%       overlap to stress-test the DoC surround-suppression mode.
%
%   Network structures (right half and bottom strip)
%     - Branching chains of overlapping capsules forming thread-like
%       networks analogous to fused/elongated mitochondria.
%
%   Rendering pipeline:
%     1. True capsule (stadium) profiles — flat body, rounded caps
%     2. Gaussian PSF blur (sigma 0.7 px) mimicking confocal optics
%     3. Poisson shot noise (120 peak photons) + Gaussian readout noise
%
% DIMENSIONS (matched to real data)
%   Typical width  : 8 px  (range 7-9 px)
%   Typical length : 12-40 px
%
% OUTPUT
%   I  – single-precision image in [0,1], size [sz x sz]

if nargin < 1, sz   = 512; end
if nargin < 2, seed = 42;  end

rng(seed);
canvas = zeros(sz, sz, 'single');

% Nominal mitochondrial dimensions
nomW = 8;    % typical width (px)
minL = 12;   % shortest rod  (px)

% -----------------------------------------------------------------------
% helpers
% -----------------------------------------------------------------------
    function canvas = addCapsule(canvas, cx, cy, len, width, angleDeg, intensity)
        % True capsule (stadium) profile:
        %   flat body along long axis, Gaussian cross-section, rounded caps
        r  = ceil(max(len/2, width/2)) + 8;
        x1 = max(1, round(cx)-r);  x2 = min(sz, round(cx)+r);
        y1 = max(1, round(cy)-r);  y2 = min(sz, round(cy)+r);
        [gx, gy] = meshgrid(single(x1:x2), single(y1:y2));
        theta = single(-deg2rad(angleDeg));
        xr    = (gx - cx).*cos(theta) - (gy - cy).*sin(theta);
        yr    = (gx - cx).*sin(theta) + (gy - cy).*cos(theta);
        halfLen      = single(len / 2);
        sigW         = single(max(width / 2.355, 0.5));
        crossProfile = exp(-(yr.^2) ./ (2 * sigW^2));
        beyond       = abs(xr) - halfLen;
        longProfile  = exp(-(max(beyond, 0).^2) ./ (2 * sigW^2));
        blob = crossProfile .* longProfile;
        canvas(y1:y2, x1:x2) = canvas(y1:y2, x1:x2) + intensity .* blob;
    end

    function canvas = addPunctum(canvas, cx, cy, radius, intensity)
        canvas = addCapsule(canvas, cx, cy, 0, radius*2, 0, intensity);
    end

% -----------------------------------------------------------------------
% 1. ISOLATED region  [1..240, 1..240]
% -----------------------------------------------------------------------

% Puncta
for k = 1:16
    cx = 20 + rand*200;  cy = 20 + rand*200;
    r  = 3 + rand*2;
    canvas = addPunctum(canvas, cx, cy, r, 0.6 + rand*0.4);
end

% Short rods — length 12-20 px
for k = 1:18
    cx  = 20 + rand*200;  cy = 20 + rand*200;
    L   = minL + rand*8;
    w   = max(nomW + randn*0.8, 6);
    ang = rand*180;
    canvas = addCapsule(canvas, cx, cy, L, w, ang, 0.6 + rand*0.4);
end

% Long rods — length 20-40 px
for k = 1:14
    cx  = 20 + rand*200;  cy = 20 + rand*200;
    L   = 20 + rand*20;
    w   = max(nomW + randn*0.8, 6);
    ang = rand*180;
    canvas = addCapsule(canvas, cx, cy, L, w, ang, 0.65 + rand*0.35);
end

% -----------------------------------------------------------------------
% 2. CLUSTERED region  [180..340, 180..340]
% -----------------------------------------------------------------------

% Cluster A — mixed puncta and short rods
centres_A = [220+randn(14,1)*16,  220+randn(14,1)*16];
for k = 1:size(centres_A,1)
    if rand < 0.4
        canvas = addPunctum(canvas, centres_A(k,1), centres_A(k,2), ...
                            3+rand*2, 0.7 + rand*0.3);
    else
        L = minL + rand*10;  w = max(nomW + randn*0.8, 6);
        canvas = addCapsule(canvas, centres_A(k,1), centres_A(k,2), ...
                            L, w, rand*180, 0.65 + rand*0.35);
    end
end

% Cluster B — roughly aligned elongated rods (fission cluster)
centres_B = [295+randn(12,1)*18,  265+randn(12,1)*18];
for k = 1:size(centres_B,1)
    L = 18 + rand*16;  w = max(nomW + randn*0.8, 6);
    canvas = addCapsule(canvas, centres_B(k,1), centres_B(k,2), ...
                        L, w, 35+randn*12, 0.65 + rand*0.35);
end

% Cluster C — very dense puncta (stress test for DoC)
centres_C = [252+randn(22,1)*10,  312+randn(22,1)*10];
for k = 1:size(centres_C,1)
    canvas = addPunctum(canvas, centres_C(k,1), centres_C(k,2), ...
                        3+rand*2, 0.55 + rand*0.45);
end

% -----------------------------------------------------------------------
% 3. NETWORK structures: right half and bottom strip
% -----------------------------------------------------------------------

% Network A — main trunk
nxt = [355, 75];  ang = -20;
for k = 1:20
    L = max(22 + randn*4, 14);  w = max(nomW + randn*0.7, 6);
    canvas = addCapsule(canvas, nxt(1), nxt(2), L, w, ang, 0.7+rand*0.3);
    ang = ang + randn*10;
    nxt = min(max(nxt + [cos(deg2rad(ang))*L*0.75, sin(deg2rad(ang))*L*0.75], 10), sz-10);
end

% Network A — branch 1
nxt = [400, 145];  ang = 55;
for k = 1:13
    L = max(20 + randn*4, 12);  w = max(nomW + randn*0.7, 6);
    canvas = addCapsule(canvas, nxt(1), nxt(2), L, w, ang, 0.65+rand*0.3);
    ang = ang + randn*13;
    nxt = min(max(nxt + [cos(deg2rad(ang))*L*0.72, sin(deg2rad(ang))*L*0.72], 10), sz-10);
end

% Network A — branch 2
nxt = [445, 210];  ang = 115;
for k = 1:14
    L = max(18 + randn*4, 12);  w = max(nomW + randn*0.7, 6);
    canvas = addCapsule(canvas, nxt(1), nxt(2), L, w, ang, 0.62+rand*0.33);
    ang = ang + randn*14;
    nxt = min(max(nxt + [cos(deg2rad(ang))*L*0.72, sin(deg2rad(ang))*L*0.72], 10), sz-10);
end

% Network B — bottom strip, tortuous
nxt = [50, 385];  ang = 8;
for k = 1:26
    L = max(20 + randn*5, 12);  w = max(nomW + randn*0.7, 6);
    canvas = addCapsule(canvas, nxt(1), nxt(2), L, w, ang, 0.65+rand*0.35);
    ang = ang + randn*16;
    nxt = min(max(nxt + [cos(deg2rad(ang))*L*0.65, sin(deg2rad(ang))*L*0.65], 10), sz-10);
end

% Network B — cross-link
nxt = [210, 428];  ang = -28;
for k = 1:15
    L = max(18 + randn*4, 12);  w = max(nomW + randn*0.7, 6);
    canvas = addCapsule(canvas, nxt(1), nxt(2), L, w, ang, 0.62+rand*0.3);
    ang = ang + randn*18;
    nxt = min(max(nxt + [cos(deg2rad(ang))*L*0.65, sin(deg2rad(ang))*L*0.65], 10), sz-10);
end

% -----------------------------------------------------------------------
% 4. PSF blur (confocal, sigma 0.7 px)
% -----------------------------------------------------------------------
psfSig = 0.7;
psf    = fspecial('gaussian', 2*ceil(3*psfSig)+1, psfSig);  % double kernel
canvas = imfilter(canvas, psf, 'replicate', 'conv');

% -----------------------------------------------------------------------
% 5. Clip, Poisson shot noise + Gaussian readout noise
% -----------------------------------------------------------------------
canvas   = min(canvas, 1);
nPhotons = 120;
counts   = single(poissrnd(double(canvas) .* nPhotons));
canvas   = counts ./ nPhotons;
canvas   = canvas + 0.012 .* single(randn(sz, sz));
I        = single(min(max(canvas, 0), 1));
end
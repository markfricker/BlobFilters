%% demoEnhancementComparison
% Demonstration script for all prototype enhancement algorithms.
%
% Generates a synthetic image with:
%   - Circular blobs of varying sizes
%   - Elongated ellipses
%   - Network-like lines
%   - Gaussian noise
%
% Then applies:
%   - logEnhance
%   - granuloLineEnhance
%   - gaussianCapsuleEnhance
%   - hessianBlobnessEnhance
%   - phaseSymmetryMinMomentEnhance
%   - fuseLogPhaseSymEnhance
%
% Segmentation is performed using LoG seeds + watershed.
%
% Prefiltering is intentionally NOT included (pipeline responsibility).

%clear; 
close all; clc;

%% ---------------- Synthetic Test Image ----------------

% rng(1);
% h = 256; w = 256;
% I = zeros(h,w);
% 
% % Circular blobs
% radii = [3 3 4 4 5 5 6 6];
% for r = radii
%     cx = randi([30 220]);
%     cy = randi([30 220]);
%     [xx,yy] = meshgrid(1:w,1:h);
%     mask = (xx-cx).^2 + (yy-cy).^2 <= r^2;
%     I(mask) = I(mask) + 1;
% end
% 
% % Elongated ellipses
% for i = 1:6
%     a = randi([8 12]);
%     b = randi([6 8]);
%     theta = rand()*pi;
%     cx = randi([40 220]);
%     cy = randi([40 220]);
%     [xx,yy] = meshgrid(1:w,1:h);
%     xRot =  (xx-cx)*cos(theta) + (yy-cy)*sin(theta);
%     yRot = -(xx-cx)*sin(theta) + (yy-cy)*cos(theta);
%     mask = (xRot.^2/a^2 + yRot.^2/b^2) <= 1;
%     I(mask) = I(mask) + 1;
% end
% 
% % Network lines
% for i = 1:10
%     x1 = randi([1 w]); y1 = randi([1 h]);
%     x2 = randi([1 w]); y2 = randi([1 h]);
%     lineMask = false(h,w);
%     idx = round(linspace(0,1,200));
%     xs = round(x1 + (x2-x1)*idx);
%     ys = round(y1 + (y2-y1)*idx);
%     valid = xs>0 & xs<=w & ys>0 & ys<=h;
%     ind = sub2ind([h w], ys(valid), xs(valid));
%     lineMask(ind) = true;
%     I(lineMask) = I(lineMask) + 0.5;
% end
% 
% % Blur + normalize
% I = imgaussfilt(I,2);
% I = I - min(I(:));
% I = I ./ max(I(:));
% 
% % Add noise
% I = imnoise(I,'gaussian',0,0.01);
% 
% %figure; imshow(I,[]); title('Synthetic Test Image');

%% ---------------- Compute LoG Seeds ----------------

seedParams.sigmas = [2 4 6];
seedParams.minDistance = 4;
seedParams.thresholdAbs = 0.15;

seedMask = computeLoGSeeds(I, seedParams);

%% ---------------- Enhancement Maps ----------------

% 1) LoG
R_log = logEnhance(I);

% 2) Granulo-Line
R_gran = granuloLineEnhance(I);

% 3) Gaussian Capsule
R_caps = gaussianCapsuleEnhance(I);

% 4) Hessian Blobness
R_hess = hessianBlobnessEnhance(I);

% 5) Phase Symmetry + Minimum Moment
[PSmap, MMmap] = phaseSymmetryMinMomentEnhance(I);

% 6) LoG + PhaseSym Fusion
fusionParams.w = 1.0;
R_fuse = fuseLogPhaseSymEnhance(I, fusionParams);

% 6) LoG + PhaseSym Fusion
R_fiber = fiberEnhance(I);

[phaseSym, orientation, totalEnergy, T, seeFilter] = phasesym_MDF(I,'polarity',1);

%% ---------------- Display Enhancement Maps ----------------

figure; 
tiledlayout(2,5)
nexttile; imshow(I,[]); title('Synthetic Test Image');
nexttile; imshow(R_log,[]); title('LoG Enhancement');
nexttile; imshow(R_gran,[]); title('Granulo-Line Enhancement');
nexttile; imshow(R_caps,[]); title('Gaussian Capsule Enhancement');
nexttile; imshow(R_hess,[]); title('Hessian Blobness Enhancement');
nexttile; imshow(PSmap,[]); title('Phase Symmetry Map');
nexttile; imshow(MMmap,[]); title('Minimum Moment Map');
nexttile; imshow(R_fuse,[]); title('LoG + PhaseSym Fusion');
nexttile; imshow(R_fiber,[]); title('fibermetric');
nexttile; imshow(phaseSym,[]); title('phase symmetry');

%% ---------------- Segmentation ----------------

segParams.thrFactor = 0.6;
segParams.minSize = 25;

L_log  = watershedSegmentation(I, R_log,  seedMask, segParams);
L_gran = watershedSegmentation(I, R_gran, seedMask, segParams);
L_caps = watershedSegmentation(I, R_caps, seedMask, segParams);
L_hess = watershedSegmentation(I, R_hess, seedMask, segParams);
L_ps   = watershedSegmentation(I, PSmap,  seedMask, segParams);
L_mm   = watershedSegmentation(I, MMmap,  seedMask, segParams);
L_fuse = watershedSegmentation(I, R_fuse, seedMask, segParams);
L_fiber = watershedSegmentation(I, R_fiber, seedMask, segParams);

%% ---------------- Display Segmentation Results ----------------

figure; 
tiledlayout(2,5)
nexttile; imshow(I,[]); title('Synthetic Test Image');
nexttile; imshow(label2rgb(L_log));  title('Segmentation: LoG');
nexttile; imshow(label2rgb(L_gran)); title('Segmentation: Granulo-Line');
nexttile; imshow(label2rgb(L_caps)); title('Segmentation: Gaussian Capsule');
nexttile; imshow(label2rgb(L_hess)); title('Segmentation: Hessian Blobness');
nexttile; imshow(label2rgb(L_ps));   title('Segmentation: Phase Symmetry');
nexttile; imshow(label2rgb(L_mm));   title('Segmentation: Minimum Moment');
nexttile; imshow(label2rgb(L_fuse)); title('Segmentation: LoG + PhaseSym');
nexttile; imshow(label2rgb(L_fiber)); title('Segmentation: fiber');

disp('Demo complete.');
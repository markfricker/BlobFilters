% demoMitoEnhance.m
%
% Demo script: synthetic + real confocal mitochondria images → five enhancers
%
%   logEnhance              – isotropic LoG bank (blob / puncta detection)
%   fiberEnhance            – fibermetric-based tubular enhancer (IPT R2018b+)
%   capsuleEnhance          – oriented capsule bank, single and DoC modes
%   granulometryEnhance     – morphological pattern spectrum (scale-integrated opening residues)
%   structureTensorEnhance  – local coherence index C ∈ [0,1]; applied as
%                             a multiplicative weight after other enhancers
%                             to separate rods (C→1) from puncta (C→0)
%
% All enhancer functions and makeMitoTestImage.m must be on the MATLAB path.
% The real image is loaded from mitImage.mat (variable I, uint16).
%
% Assumed mitochondrial dimensions (pixels):
%   width  : ~8 px
%   length : 12-40 px

clear; clc; close all;

% =========================================================================
% 1.  Load / generate images
% =========================================================================
fprintf('Generating synthetic test image...\n');
Isynth = makeMitoTestImage(512, 42);

fprintf('Loading real image...\n');
realData = load('mitImage.mat');
Ireal    = im2single(realData.I);   % uint16 → single [0,1]

% =========================================================================
% 2.  Enhancer parameters  (tuned to width~8px, length 12-40px)
% =========================================================================

% --- logEnhance -----------------------------------------------------------
pLog.sigmas    = [2 3 4 5 6];      % sigma~4 → peak response at dia~8px
pLog.normalize = true;

% --- fiberEnhance ---------------------------------------------------------
pFib.widths    = [6 7 8 9 10];     % centred on nominal 8px width
pFib.multimode = 'stack';
pFib.normalize = true;

% --- capsuleEnhance : single mode ----------------------------------------
pCapS.lengths      = [12 16 20 28 36 40];
pCapS.width        = 8;
pCapS.orientations = 12;
pCapS.mode         = 'single';
pCapS.normalize    = true;

% --- capsuleEnhance : DoC mode -------------------------------------------
pCapD.lengths      = [12 16 20 28 36 40];
pCapD.width        = 8;
pCapD.wideWidth    = 18;           % inhibitory surround ~2.25x nominal width
pCapD.alpha        = 0.55;
pCapD.orientations = 12;
pCapD.mode         = 'doc';
pCapD.normalize    = true;

% --- granulometryEnhance --------------------------------------------------
% Radii are chosen to span the mitochondrial cross-section (half-width to
% full-width) and capture both puncta (~4 px radius) and wider tubules.
% The pattern spectrum differentiates between successive disk openings, so
% radii should be roughly evenly spaced in the range of interest.
% Upper bound (r=16) is set to suppress larger background structures.
pGran.sigmas    = [2 4 6 8 10 12 16];   % disk radii in px; width~8px → peak at r=4-8
pGran.normalize = true;

% --- structureTensorEnhance -----------------------------------------------
% sigmaGrad: gradient smoothing scale — suppress pixel-level noise.
% sigmaInt : integration scale — should match structure half-width (~4-5 px).
% normalize: rescale C to [0,1] so it acts as a clean [0,1] weight.
pST.sigmaGrad  = 1.5;
pST.sigmaInt   = 5;
pST.normalize  = true;

% =========================================================================
% 3.  Run enhancers on BOTH images
% =========================================================================
images  = {Isynth, 'Synthetic'; Ireal, 'Real (mitImage)'};
results = cell(2, 6);   % {logR, fibR, capSR, capDR, granR, coherR} per image

for im = 1:2
    img = images{im,1};
    lbl = images{im,2};
    fprintf('\n--- %s ---\n', lbl);

    fprintf('  logEnhance...         '); tic;
    results{im,1} = logEnhance(img, pLog);
    fprintf('%.2fs\n', toc);

    fprintf('  fiberEnhance...       ');
    try
        tic; results{im,2} = fiberEnhance(img, pFib); fprintf('%.2fs\n', toc);
    catch ME
        fprintf('SKIPPED (%s)\n', ME.message);
        results{im,2} = zeros(size(img), 'single');
    end

    fprintf('  capsule single...     '); tic;
    results{im,3} = capsuleEnhance(img, pCapS);
    fprintf('%.2fs\n', toc);

    fprintf('  capsule DoC...        '); tic;
    results{im,4} = capsuleEnhance(img, pCapD);
    fprintf('%.2fs\n', toc);

    fprintf('  granulometryEnhance.. '); tic;
    results{im,5} = granulometryEnhance(img, pGran);
    fprintf('%.2fs\n', toc);

    fprintf('  structureTensor...    '); tic;
    results{im,6} = structureTensorEnhance(img, pST);
    fprintf('%.2fs\n', toc);
end

% =========================================================================
% 4.  Figure 1 — annotated synthetic test image
% =========================================================================
figure(1);
set(gcf,'Name','Synthetic Test Image','NumberTitle','off', ...
        'Color','k','Position',[30 30 560 560]);
imshow(Isynth, []);  colormap(gca,'gray');
title('Synthetic Mitochondria Test Image','Color','w','FontSize',13);
hold on;
rSpec = {[10  10  230 230],'Isolated',  [0.4 0.8 0.4];
         [175 175 175 175],'Clustered', [0.9 0.7 0.2];
         [330  10  170 330],'Network A', [0.3 0.6 0.9];
         [10  360  490 140],'Network B', [0.9 0.4 0.4]};
for k = 1:size(rSpec,1)
    b = rSpec{k,1};
    rectangle('Position',b,'EdgeColor',rSpec{k,3},'LineWidth',1.5,'LineStyle','--');
    text(b(1)+4, b(2)+14, rSpec{k,2}, 'Color',rSpec{k,3}, ...
         'FontSize',9,'FontWeight','bold');
end
hold off;

% =========================================================================
% 5.  Figure 2 — synthetic image comparison panel (now 6 panels, 2x3)
% =========================================================================
titles2 = {'Raw (synthetic)', ...
           'logEnhance', ...
           'fiberEnhance', ...
           'capsule single', ...
           'capsule DoC', ...
           'granulometry'};
panels2 = {Isynth, results{1,1}, results{1,2}, results{1,3}, results{1,4}, results{1,5}};

figure(2);
set(gcf,'Name','Synthetic — Enhancement Comparison','NumberTitle','off', ...
        'Color','k','Position',[30 30 1250 840]);
cmaps = {'gray','hot','hot','hot','hot','hot'};
for k = 1:6
    ax = subplot(2,3,k);
    imshow(panels2{k},[]); colormap(ax, cmaps{k});
    title(titles2{k},'Color','w','FontSize',10);
    set(ax,'XColor','none','YColor','none');
end
sgtitle('Synthetic Image — Enhancement Comparison', ...
        'Color','w','FontSize',13,'FontWeight','bold');

% =========================================================================
% 6.  Figure 3 — real image comparison panel (now 6 panels, 2x3)
% =========================================================================
titles3 = {'Raw (real)', ...
           'logEnhance', ...
           'fiberEnhance', ...
           'capsule single', ...
           'capsule DoC', ...
           'granulometry'};
panels3 = {Ireal, results{2,1}, results{2,2}, results{2,3}, results{2,4}, results{2,5}};

figure(3);
set(gcf,'Name','Real Image — Enhancement Comparison','NumberTitle','off', ...
        'Color','k','Position',[60 60 1250 840]);
for k = 1:6
    ax = subplot(2,3,k);
    imshow(panels3{k},[]); colormap(ax, cmaps{k});
    title(titles3{k},'Color','w','FontSize',10);
    set(ax,'XColor','none','YColor','none');
end
sgtitle('Real Image — Enhancement Comparison', ...
        'Color','w','FontSize',13,'FontWeight','bold');

% =========================================================================
% 7.  Figure 4 — DoC vs single difference map + cluster close-up (synthetic)
% =========================================================================
diff_DS = results{1,4} - results{1,3};

figure(4);
set(gcf,'Name','DoC vs Single','NumberTitle','off', ...
        'Color','k','Position',[90 90 1100 460]);

% Build blue-white-red diverging colormap
n   = 256;  t = linspace(0,1,n)';
bwr = [t, t, ones(n,1); ones(n,1), flipud(t), flipud(t)];
bwr = bwr(round(linspace(1,size(bwr,1),n)),:);

ax1 = subplot(1,3,1);
imshow(diff_DS,[-0.5 0.5]); colormap(ax1,bwr);
cb = colorbar; cb.Color = 'w';
title('DoC − Single (synthetic)','Color','w','FontSize',10);
set(ax1,'XColor','none','YColor','none');

roi = [180 180 160 160];
ax2 = subplot(1,3,2);
imshow(imcrop(results{1,3},roi),[]); colormap(ax2,'hot');
title('Cluster close-up: single','Color','w','FontSize',10);
set(ax2,'XColor','none','YColor','none');

ax3 = subplot(1,3,3);
imshow(imcrop(results{1,4},roi),[]); colormap(ax3,'hot');
title('Cluster close-up: DoC','Color','w','FontSize',10);
set(ax3,'XColor','none','YColor','none');

sgtitle('DoC surround-suppression on clustered region', ...
        'Color','w','FontSize',13,'FontWeight','bold');

% =========================================================================
% 8.  Figure 5 — DoC vs single on real image
% =========================================================================
figure(5);
set(gcf,'Name','DoC vs Single — Real','NumberTitle','off', ...
        'Color','k','Position',[120 120 900 420]);

ax4 = subplot(1,2,1);
imshow(results{2,3},[]); colormap(ax4,'hot');
title('Real: capsule single','Color','w','FontSize',11);
set(ax4,'XColor','none','YColor','none');

ax5 = subplot(1,2,2);
imshow(results{2,4},[]); colormap(ax5,'hot');
title('Real: capsule DoC','Color','w','FontSize',11);
set(ax5,'XColor','none','YColor','none');

sgtitle('Real Image — Single vs DoC', ...
        'Color','w','FontSize',13,'FontWeight','bold');

% =========================================================================
% 9.  Figure 6 — granulometry vs logEnhance comparison
%     Both methods integrate responses across spatial scales; this figure
%     shows where the morphological (granulometry) and Laplacian-of-Gaussian
%     approaches agree and diverge on the same scene.
% =========================================================================
diff_gran_log_synth = results{1,5} - results{1,1};
diff_gran_log_real  = results{2,5} - results{2,1};

figure(6);
set(gcf,'Name','Granulometry vs LoG','NumberTitle','off', ...
        'Color','k','Position',[150 150 1250 840]);

% Row 1: synthetic
ax6 = subplot(2,4,1);
imshow(Isynth,[]); colormap(ax6,'gray');
title('Raw (synthetic)','Color','w','FontSize',10);
set(ax6,'XColor','none','YColor','none');

ax7 = subplot(2,4,2);
imshow(results{1,1},[]); colormap(ax7,'hot');
title('logEnhance','Color','w','FontSize',10);
set(ax7,'XColor','none','YColor','none');

ax8 = subplot(2,4,3);
imshow(results{1,5},[]); colormap(ax8,'hot');
title('granulometry','Color','w','FontSize',10);
set(ax8,'XColor','none','YColor','none');

ax9 = subplot(2,4,4);
imshow(diff_gran_log_synth,[-0.5 0.5]); colormap(ax9,bwr);
cb2 = colorbar; cb2.Color = 'w';
title('Gran − LoG (synthetic)','Color','w','FontSize',10);
set(ax9,'XColor','none','YColor','none');

% Row 2: real
ax10 = subplot(2,4,5);
imshow(Ireal,[]); colormap(ax10,'gray');
title('Raw (real)','Color','w','FontSize',10);
set(ax10,'XColor','none','YColor','none');

ax11 = subplot(2,4,6);
imshow(results{2,1},[]); colormap(ax11,'hot');
title('logEnhance','Color','w','FontSize',10);
set(ax11,'XColor','none','YColor','none');

ax12 = subplot(2,4,7);
imshow(results{2,5},[]); colormap(ax12,'hot');
title('granulometry','Color','w','FontSize',10);
set(ax12,'XColor','none','YColor','none');

ax13 = subplot(2,4,8);
imshow(diff_gran_log_real,[-0.5 0.5]); colormap(ax13,bwr);
cb3 = colorbar; cb3.Color = 'w';
title('Gran − LoG (real)','Color','w','FontSize',10);
set(ax13,'XColor','none','YColor','none');

sgtitle('Granulometry vs LoG — Scale-integrated morphological vs Laplacian response', ...
        'Color','w','FontSize',13,'FontWeight','bold');

% =========================================================================
% 10. Figure 7 — Structure tensor coherence: rod vs puncta separation
%
%     Coherence C ≈ 1 for elongated rods/fibres, C ≈ 0 for isotropic
%     puncta and noise.  Applied as a multiplicative weight:
%
%       DoC × C       → rod / tubular channel   (puncta suppressed)
%       LoG × (1−C)   → puncta channel          (rods suppressed)
%
%     Both derived maps are individually rescaled to [0,1] for display.
% =========================================================================
figure(7);
set(gcf,'Name','Structure Tensor Coherence','NumberTitle','off', ...
        'Color','k','Position',[180 180 1250 840]);

for im = 1:2
    C      = results{im,6};          % coherence map
    capDoc = results{im,4};          % capsule DoC — rod-selective enhancer
    logR   = results{im,1};          % logEnhance  — puncta-inclusive enhancer
    lbl    = images{im,2};

    rodW    = capDoc .* C;
    mx = max(rodW(:));   if mx > 0, rodW    = rodW    / mx; end

    punctaW = logR .* (1 - C);
    mx = max(punctaW(:)); if mx > 0, punctaW = punctaW / mx; end

    base = (im-1)*4;

    ax = subplot(2,4, base+1);
    imshow(images{im,1},[]); colormap(ax,'gray');
    title(sprintf('Raw (%s)', lbl),'Color','w','FontSize',9);
    set(ax,'XColor','none','YColor','none');

    ax = subplot(2,4, base+2);
    imshow(C,[]); colormap(ax,'parula');
    title('Coherence C','Color','w','FontSize',9);
    set(ax,'XColor','none','YColor','none');

    ax = subplot(2,4, base+3);
    imshow(rodW,[]); colormap(ax,'hot');
    title('DoC \times C  (rods)','Color','w','FontSize',9);
    set(ax,'XColor','none','YColor','none');

    ax = subplot(2,4, base+4);
    imshow(punctaW,[]); colormap(ax,'hot');
    title('LoG \times (1-C)  (puncta)','Color','w','FontSize',9);
    set(ax,'XColor','none','YColor','none');
end
sgtitle('Structure Tensor Coherence — Rod vs Puncta separation', ...
        'Color','w','FontSize',13,'FontWeight','bold');

fprintf('\nDone. Seven figures generated.\n');
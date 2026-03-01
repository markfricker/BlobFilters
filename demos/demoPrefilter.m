% demoPrefilter.m
%
% Demo script: compare three pre-filtering strategies and their effect on
% downstream rod enhancement (capsule DoC) and structure tensor coherence.
%
%   imdiffusefilt     – Perona-Malik anisotropic diffusion (IPT built-in)
%                       Isotropic edge-stopping: smooths within regions,
%                       stops at intensity edges.  Not directionally steered.
%
%   imguidedfilter    – Guided filter (IPT built-in, R2014a+)
%                       Intensity-similarity gated edge-preserving smoothing.
%                       Fast (linear complexity); not directionally steered.
%
%   orientedGaussSmooth – Orientation-adaptive Gaussian smoothing (this toolbox)
%                       Steers the smoothing kernel along the local fibre
%                       direction using the structure tensor.  One-pass
%                       non-iterative approximation to CED.
%
% WHAT THE FIGURES SHOW
%   For each image (synthetic, real) a 3×4 panel is produced:
%     Row 1 – pre-filtered image under each method (visual quality)
%     Row 2 – capsule DoC response applied to each pre-filtered image
%     Row 3 – structure tensor coherence C from each pre-filtered image
%
%   A cleaner C map (fewer spurious high-C pixels in noise/background)
%   indicates that the pre-filter improved orientation estimation.
%   Brighter DoC responses on true rods with fewer false responses on
%   puncta indicate better rod selectivity.
%
% PARAMETERS (tuned to width~8px, length 12-40px)
%   imdiffusefilt:     NumberOfIterations=15, GradientThreshold=0.03, ConductionMethod='exponential'
%   imguidedfilter:    NeighborhoodSize=7, DegreeOfSmoothing=5e-3
%   orientedGaussSmooth: sigmaAlong=4, sigmaAcross=1.5, orientations=8

clear; clc; close all;

% =========================================================================
% 1.  Load / generate images
% =========================================================================
fprintf('Generating synthetic test image...\n');
Isynth = makeMitoTestImage(512, 42);

fprintf('Loading real image...\n');
realData = load('mitImage.mat');
Ireal    = im2single(realData.I);

% =========================================================================
% 2.  Pre-filter parameters
% =========================================================================

% --- Perona-Malik anisotropic diffusion -----------------------------------
% GradientThreshold: gradients above this value stop diffusion.
% For [0,1] single images, mitochondria boundary gradients ~0.02-0.1;
% 0.03 preserves boundaries while smoothing noise.
pmIter  = 15;
pmGrad  = 0.03;

% --- Guided filter --------------------------------------------------------
% NeighborhoodSize: local window (px); should be ~half the structure width.
% DegreeOfSmoothing: ε² regularisation; smaller = more edge-preserving.
gfNbhd  = 7;
gfEps   = 5e-3;

% --- Orientation-adaptive Gaussian smoothing ------------------------------
pOGS.sigmaAlong   = 4;     % smooth 4 px along fibres
pOGS.sigmaAcross  = 1.5;   % smooth 1.5 px across fibres (preserve edges)
pOGS.orientations = 8;
pOGS.sigmaGrad    = 1.5;
pOGS.sigmaInt     = 5;

% =========================================================================
% 3.  Downstream enhancer parameters (capsule DoC + structure tensor)
% =========================================================================
pCapD.lengths      = [12 16 20 28 36 40];
pCapD.width        = 8;
pCapD.wideWidth    = 18;
pCapD.alpha        = 0.55;
pCapD.orientations = 12;
pCapD.mode         = 'doc';
pCapD.normalize    = true;

pST.sigmaGrad  = 1.5;
pST.sigmaInt   = 5;
pST.normalize  = true;

% =========================================================================
% 4.  Apply pre-filters and run downstream on BOTH images
% =========================================================================
% Storage: {raw, pm, guided, ogs} × {filtered_img, capDoc, coherence}
images = {Isynth, 'Synthetic'; Ireal, 'Real (mitImage)'};

for im = 1:2
    img = images{im,1};
    lbl = images{im,2};
    fprintf('\n--- %s ---\n', lbl);

    % ---- pre-filters -------------------------------------------------------
    fprintf('  imdiffusefilt...      '); tic;
    I_pm  = imdiffusefilt(img, ...
                'NumberOfIterations', pmIter, ...
                'GradientThreshold',  pmGrad, ...
                'ConductionMethod',   'exponential');
    fprintf('%.2fs\n', toc);

    fprintf('  imguidedfilter...     '); tic;
    I_gf  = imguidedfilter(img, img, ...
                'NeighborhoodSize',   gfNbhd, ...
                'DegreeOfSmoothing',  gfEps);
    fprintf('%.2fs\n', toc);

    fprintf('  orientedGaussSmooth.. '); tic;
    I_ogs = orientedGaussSmooth(img, pOGS);
    fprintf('%.2fs\n', toc);

    filtered = {img, I_pm, I_gf, I_ogs};   % raw + 3 filtered

    % ---- downstream enhancers on each filtered image -----------------------
    capD    = cell(1, 4);
    coher   = cell(1, 4);
    labels  = {'Raw', 'Perona-Malik', 'Guided', 'Oriented Gauss'};

    for f = 1:4
        fi = filtered{f};
        fprintf('  capsule DoC [%s]...', labels{f}); tic;
        capD{f}  = capsuleEnhance(fi, pCapD);
        fprintf('%.2fs\n', toc);

        coher{f} = structureTensorEnhance(fi, pST);
    end

    % =========================================================================
    % 5.  Figure (2 per image): 3 rows × 4 cols
    %     Row 1 – pre-filtered images
    %     Row 2 – capsule DoC response
    %     Row 3 – structure tensor coherence C
    % =========================================================================
    figNum = im;   % Fig 1 = synthetic, Fig 2 = real
    figure(figNum);
    set(gcf,'Name', sprintf('Pre-filter comparison — %s', lbl), ...
            'NumberTitle','off','Color','k', ...
            'Position', [30 + (im-1)*40, 30 + (im-1)*40, 1400, 900]);

    rowLabel = {'Pre-filtered', 'Capsule DoC', 'Coherence C'};
    cmaps    = {{'gray','gray','gray','gray'}, ...
                {'hot','hot','hot','hot'}, ...
                {'parula','parula','parula','parula'}};

    for row = 1:3
        switch row
            case 1, data = filtered;
            case 2, data = capD;
            case 3, data = coher;
        end

        for col = 1:4
            ax = subplot(3, 4, (row-1)*4 + col);
            imshow(data{col}, []); colormap(ax, cmaps{row}{col});
            if row == 1
                title(labels{col}, 'Color','w','FontSize',10,'FontWeight','bold');
            end
            if col == 1
                ylabel(rowLabel{row}, 'Color','w','FontSize',9);
            end
            set(ax,'XColor','none','YColor','none');
        end
    end

    sgtitle(sprintf('Pre-filter comparison — %s image', lbl), ...
            'Color','w','FontSize',13,'FontWeight','bold');
end

fprintf('\nDone. Two figures generated.\n');

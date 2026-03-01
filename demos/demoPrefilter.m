% demoPrefilter.m
%
% Demo script: compare three pre-filtering strategies and their effect on
% the two most powerful downstream rod enhancers.
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
%     Row 2 – rodGranulometryEnhance response on each pre-filtered image
%     Row 3 – fiberEnhance response on each pre-filtered image
%
%   A brighter response on true rods/fibres with less background and fewer
%   false responses on puncta indicates the pre-filter improved detection.
%   rodGranulometryEnhance (morphological, line SEs) and fiberEnhance
%   (Hessian-based) are complementary: comparing both shows whether the
%   improvement is consistent across detector types.
%
% PARAMETERS (tuned to width~8px, length 12-40px)
%   imdiffusefilt:       NumberOfIterations=15, GradientThreshold=0.03, ConductionMethod='exponential'
%   imguidedfilter:      NeighborhoodSize=7, DegreeOfSmoothing=5e-3
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
% 3.  Downstream enhancer parameters (rodGranulometry + fiberEnhance)
% =========================================================================

% --- rodGranulometryEnhance -----------------------------------------------
pRodGran.lengths      = [8 12 16 20 28 36];
pRodGran.orientations = 8;
pRodGran.normalize    = true;

% --- fiberEnhance ---------------------------------------------------------
% fibermetric is available in IPT R2018b+; the try/catch below handles
% older toolbox versions gracefully.
pFib.widths    = [6 7 8 9 10];
pFib.multimode = 'stack';
pFib.normalize = true;

% =========================================================================
% 4.  Apply pre-filters and run downstream on BOTH images
% =========================================================================
% Storage: {raw, pm, guided, ogs} × {filtered_img, rodGran, fiberR}
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
    labels   = {'Raw', 'Perona-Malik', 'Guided', 'Oriented Gauss'};

    % ---- downstream enhancers on each filtered image -----------------------
    rodGranR = cell(1, 4);
    fiberR   = cell(1, 4);

    for f = 1:4
        fi = filtered{f};

        fprintf('  rodGranulometry [%s]... ', labels{f}); tic;
        rodGranR{f} = rodGranulometryEnhance(fi, pRodGran);
        fprintf('%.2fs\n', toc);

        fprintf('  fiberEnhance    [%s]... ', labels{f});
        try
            tic; fiberR{f} = fiberEnhance(fi, pFib); fprintf('%.2fs\n', toc);
        catch ME
            fprintf('SKIPPED (%s)\n', ME.message);
            fiberR{f} = zeros(size(fi), 'single');
        end
    end

    % =========================================================================
    % 5.  Figure (2 per image): 3 rows × 4 cols
    %     Row 1 – pre-filtered images
    %     Row 2 – rodGranulometryEnhance response
    %     Row 3 – fiberEnhance response
    % =========================================================================
    figNum = im;   % Fig 1 = synthetic, Fig 2 = real
    figure(figNum);
    set(gcf,'Name', sprintf('Pre-filter comparison — %s', lbl), ...
            'NumberTitle','off','Color','k', ...
            'Position', [30 + (im-1)*40, 30 + (im-1)*40, 1400, 900]);

    rowLabel = {'Pre-filtered', 'Rod granulometry', 'Fiber enhance'};
    cmaps    = {{'gray','gray','gray','gray'}, ...
                {'hot','hot','hot','hot'}, ...
                {'hot','hot','hot','hot'}};

    for row = 1:3
        switch row
            case 1, data = filtered;
            case 2, data = rodGranR;
            case 3, data = fiberR;
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

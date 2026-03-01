% demoPrefilter.m
%
% Demo script: compare three pre-filtering strategies and their effect on
% the two most powerful downstream rod enhancers.
%
%   imdiffusefilt       – Perona-Malik anisotropic diffusion (IPT built-in)
%                         Isotropic edge-stopping: smooths within regions,
%                         stops at intensity edges.  Not directionally steered.
%
%   imguidedfilter      – Guided filter (IPT built-in, R2014a+)
%                         Intensity-similarity gated edge-preserving smoothing.
%                         Fast (linear complexity); not directionally steered.
%
%   orientedGaussSmooth – Orientation-adaptive Gaussian smoothing (this toolbox)
%                         Steers the smoothing kernel along the local fibre
%                         direction using the structure tensor.  One-pass
%                         non-iterative approximation to CED.
%
% WHAT THE FIGURES SHOW
%   One figure per pre-filter method (4 total, including Raw baseline), each
%   a 2×3 panel:
%     Rows – Synthetic image (top) / Real image (bottom)
%     Cols – Pre-filtered image | rodGranulometryEnhance | fiberEnhance
%
%   Comparing across figures shows how much each pre-filter improves
%   detection; rodGranulometryEnhance (morphological) and fiberEnhance
%   (Hessian-based) are complementary detectors so consistent improvement
%   in both columns is a stronger indicator.
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

images   = {Isynth, 'Synthetic'; Ireal, 'Real (mitImage)'};

% =========================================================================
% 2.  Pre-filter parameters
% =========================================================================

% --- Perona-Malik anisotropic diffusion -----------------------------------
pmIter = 15;
pmGrad = 0.03;

% --- Guided filter --------------------------------------------------------
gfNbhd = 7;
gfEps  = 5e-3;

% --- Orientation-adaptive Gaussian smoothing ------------------------------
pOGS.sigmaAlong   = 4;
pOGS.sigmaAcross  = 1.5;
pOGS.orientations = 8;
pOGS.sigmaGrad    = 1.5;
pOGS.sigmaInt     = 5;

% =========================================================================
% 3.  Downstream enhancer parameters
% =========================================================================

% --- rodGranulometryEnhance -----------------------------------------------
pRodGran.lengths      = [8 12 16 20 28 36];
pRodGran.orientations = 8;
pRodGran.normalize    = true;

% --- fiberEnhance (IPT R2018b+) ------------------------------------------
pFib.widths    = [6 7 8 9 10];
pFib.multimode = 'stack';
pFib.normalize = true;

% =========================================================================
% 4.  Pre-filter all images
%     filtered{im, method}  im=1 synthetic, im=2 real
%                           method=1 raw, 2 PM, 3 guided, 4 OGS
% =========================================================================
methodNames = {'Raw (no pre-filter)', 'Perona-Malik', 'Guided', 'Oriented Gauss'};
nMethods    = numel(methodNames);

filtered = cell(2, nMethods);

for im = 1:2
    img = images{im,1};
    lbl = images{im,2};
    fprintf('\n--- Pre-filtering: %s ---\n', lbl);

    filtered{im,1} = img;   % raw baseline

    fprintf('  imdiffusefilt...      '); tic;
    filtered{im,2} = imdiffusefilt(img, ...
        'NumberOfIterations', pmIter, ...
        'GradientThreshold',  pmGrad, ...
        'ConductionMethod',   'exponential');
    fprintf('%.2fs\n', toc);

    fprintf('  imguidedfilter...     '); tic;
    filtered{im,3} = imguidedfilter(img, img, ...
        'NeighborhoodSize',  gfNbhd, ...
        'DegreeOfSmoothing', gfEps);
    fprintf('%.2fs\n', toc);

    fprintf('  orientedGaussSmooth.. '); tic;
    filtered{im,4} = orientedGaussSmooth(img, pOGS);
    fprintf('%.2fs\n', toc);
end

% =========================================================================
% 5.  Run downstream enhancers on all filtered images
%     rodGranR{im, method}, fiberR{im, method}
% =========================================================================
rodGranR = cell(2, nMethods);
fiberR   = cell(2, nMethods);

for im = 1:2
    lbl = images{im,2};
    fprintf('\n--- Downstream enhancers: %s ---\n', lbl);

    for f = 1:nMethods
        fi = filtered{im,f};

        fprintf('  rodGranulometry [%-14s]  ', methodNames{f}); tic;
        rodGranR{im,f} = rodGranulometryEnhance(fi, pRodGran);
        fprintf('%.2fs\n', toc);

        fprintf('  fiberEnhance    [%-14s]  ', methodNames{f});
        try
            tic; fiberR{im,f} = fiberEnhance(fi, pFib); fprintf('%.2fs\n', toc);
        catch ME
            fprintf('SKIPPED (%s)\n', ME.message);
            fiberR{im,f} = zeros(size(fi), 'single');
        end
    end
end

% =========================================================================
% 6.  Figures — one per pre-filter method
%     Layout: 2 rows (Synthetic top, Real bottom) × 3 cols
%             Col 1: pre-filtered image (gray)
%             Col 2: rodGranulometryEnhance  (hot)
%             Col 3: fiberEnhance            (hot)
% =========================================================================
colTitles = {'Pre-filtered', 'Rod granulometry', 'Fiber enhance'};
rowLabels = {images{1,2}, images{2,2}};

for f = 1:nMethods
    figure(f);
    set(gcf, 'Name', methodNames{f}, 'NumberTitle','off', ...
             'Color','k', 'Position', [30 + (f-1)*30, 30 + (f-1)*30, 1100, 740]);

    for im = 1:2
        panelData = {filtered{im,f}, rodGranR{im,f}, fiberR{im,f}};
        cmaps     = {'gray', 'hot', 'hot'};

        for col = 1:3
            ax = subplot(2, 3, (im-1)*3 + col);
            imshow(panelData{col}, []); colormap(ax, cmaps{col});

            % Column titles on top row only
            if im == 1
                title(colTitles{col}, 'Color','w', 'FontSize',10);
            end
            % Row labels on left column only
            if col == 1
                ylabel(rowLabels{im}, 'Color','w', 'FontSize',10);
            end
            set(ax, 'XColor','none', 'YColor','none');
        end
    end

    sgtitle(sprintf('Pre-filter: %s', methodNames{f}), ...
            'Color','w', 'FontSize',13, 'FontWeight','bold');
end

fprintf('\nDone. %d figures generated (one per method).\n', nMethods);

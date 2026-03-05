%% 
% demoMitoEnhance.m
%
% Demo script: classical enhancers (OGS-smoothed + coherence-weighted)
%              versus Cellpose deep-learning segmentation
%
% PIPELINE
%   Preprocessing  : orientedGaussSmooth(Ireal)  → Ism
%                    Smooths along the local fibre axis; suppresses noise
%                    without blurring mitochondrial boundaries.
%   Classical enhancers (run on Ism):
%     logEnhance
%     fiberEnhance
%     capsuleEnhance  (DoC)
%     rodGranulometryEnhance
%   Deep-learning  : cellposeEnhance(Ireal)  → binary mask + label image
%                    Run on the original image (Cellpose normalises internally).
%
% FIGURES PRODUCED (real image only; synthetic toggled OFF)
%   Figure 1 — full field of view:
%     Raw | log | fiber | capDoC | rodGran | CP mask | CP labels
%   Figure 2 — zoom roi1 (cluster 1):  same seven columns
%   Figure 3 — zoom roi2 (cluster 2):  same seven columns
%
% REQUIREMENTS
%   All .m enhancer files + makeMitoTestImage.m on the MATLAB path
%   mitImage.mat  (variable I, uint16 confocal mito image)
%   For Cellpose: Medical Imaging Toolbox Interface for Cellpose Library add-on
%                 + Python cellpose >= 3.x in the configured pyenv environment

clear; clc; close all;

% =========================================================================
% TOGGLE
% =========================================================================
showSynthetic = false;   % set true to re-enable synthetic-image figures

% =========================================================================
% ZOOM REGIONS  ← adjust after Figure 1 to frame the two main clusters
%
%   Format: [x_left, y_top, width, height]  (imcrop convention, 1-based px)
%   Tip:    hover over Figure 1 with the data cursor to read pixel coords.
% =========================================================================
roi1 = [100  80  220 220];   % cluster 1  (adjust)
roi2 = [330 220  220 220];   % cluster 2  (adjust)

% =========================================================================
% 1.  Load images
% =========================================================================
if showSynthetic
    fprintf('Generating synthetic test image...\n');
    Isynth = makeMitoTestImage(512, 42);
end

fprintf('Loading real image...\n');
realData = load('mitImage.mat');
Ireal    = im2single(realData.I);   % uint16 → single [0,1]
[imgH, imgW] = size(Ireal);
fprintf('  Image size: %d × %d px\n', imgW, imgH);

% =========================================================================
% 2.  Parameters
% =========================================================================

% --- orientedGaussSmooth (preprocessing) ---------------------------------
pOGS.sigmaAlong   = 4;    % along-fibre smoothing scale (px)
pOGS.sigmaAcross  = 1.5;  % across-fibre smoothing scale (px)
pOGS.orientations = 8;
pOGS.sigmaGrad    = 1.5;
pOGS.sigmaInt     = 5;

% --- logEnhance ----------------------------------------------------------
pLog.sigmas    = [2 3 4 5 6];
pLog.normalize = true;

% --- fiberEnhance --------------------------------------------------------
pFib.widths    = [6 7 8 9 10];
pFib.multimode = 'stack';
pFib.normalize = true;

% --- capsuleEnhance DoC --------------------------------------------------
pCapD.lengths      = [12 16 20 28 36 40];
pCapD.width        = 8;
pCapD.wideWidth    = 18;
pCapD.alpha        = 0.55;
pCapD.orientations = 12;
pCapD.mode         = 'doc';
pCapD.normalize    = true;

% --- rodGranulometryEnhance ----------------------------------------------
pRodGran.lengths      = [8 12 16 20 28 36];
pRodGran.orientations = 8;
pRodGran.normalize    = true;

% --- cellposeEnhance -----------------------------------------------------
% cyto3: Cellpose 3 super-generalist model (Nature Methods 2025).
%   Requires Python cellpose >= 3.x; weights auto-downloaded on first call.
% cellProb = -2: relaxed threshold to recover faint / elongated mitochondria
%   (default 0; range -6 to 6 — lower = more detections, more noise below -4)
% flowThreshold = 0.8: tolerate imperfect flow fields on elongated objects
%   (default 0.4; range 0.1-3 — higher = more detections, rougher boundaries)
pCP.model         = 'cyto3';
pCP.diameter      = 10;
pCP.cellProb      = 0;
pCP.flowThreshold = 0.8;

% =========================================================================
% 3.  Preprocessing: OGS smooth
% =========================================================================
fprintf('\n--- Preprocessing ---\n');

fprintf('  orientedGaussSmooth...    '); tic;
Ism = orientedGaussSmooth(Ireal, pOGS);
fprintf('%.2fs\n', toc);

% =========================================================================
% 4.  Classical enhancers on Ism
% =========================================================================
fprintf('\n--- Classical enhancers (on OGS-smoothed image) ---\n');

fprintf('  logEnhance...             '); tic;
logR  = logEnhance(Ism, pLog);
fprintf('%.2fs\n', toc);

fprintf('  fiberEnhance...           ');
try
    tic; fibR = fiberEnhance(Ism, pFib); fprintf('%.2fs\n', toc);
catch ME
    fprintf('SKIPPED (%s)\n', ME.message);
    fibR = zeros(size(Ism), 'single');
end

fprintf('  capsule DoC...            '); tic;
capDR = capsuleEnhance(Ism, pCapD);
fprintf('%.2fs\n', toc);

fprintf('  rodGranulometry...        '); tic;
rodR  = rodGranulometryEnhance(Ism, pRodGran);
fprintf('%.2fs\n', toc);

% (enhancers already output normalised [0,1] maps via normalize=true)

% =========================================================================
% 5.  Cellpose on original image
% =========================================================================
fprintf('\n--- Cellpose ---\n');
fprintf('  cellposeEnhance...        ');
hasCellpose = exist('cellpose', 'file') ~= 0;
if hasCellpose
    try
        tic;
        [cpMask, cpL] = cellposeEnhance(Ireal, pCP);
        fprintf('%.2fs  (%d objects)\n', toc, max(cpL(:)));
    catch ME
        fprintf('SKIPPED (%s)\n', ME.message);
        cpMask      = zeros(size(Ireal), 'single');
        cpL         = zeros(size(Ireal), 'uint16');
        hasCellpose = false;
    end
else
    fprintf('SKIPPED (Cellpose add-on not installed)\n');
    cpMask = zeros(size(Ireal), 'single');
    cpL    = zeros(size(Ireal), 'uint16');
end

% Coloured instance-label overlay (uint8 RGB)
if hasCellpose && max(cpL(:)) > 0
    cpRGB = label2rgb(cpL, 'hsv', 'k');
else
    cpRGB = zeros(imgH, imgW, 3, 'uint8');
end

% =========================================================================
% 6.  Panel layout (seven columns, constant across all three figures)
% =========================================================================
nLabel      = max(cpL(:));
panelTitles = {'Raw', ...
               'log', ...
               'fiber', ...
               'capDoC', ...
               'rodGran', ...
               sprintf('Cellpose mask  (n=%d)', nLabel), ...
               sprintf('Cellpose labels  (n=%d)', nLabel)};
cmaps = {'gray', 'hot', 'hot', 'hot', 'hot', 'hot', []};
% [] = truecolor RGB for labels panel

% =========================================================================
% 7.  Figure 1 — Full FOV
% =========================================================================
figure(1);
set(gcf, 'Name', 'Real Image — Full FOV', 'NumberTitle', 'off', ...
         'Color', 'k', 'Position', [30 680 1820 380]);

mitoFillFigure(1, ...
    {Ireal, logR, fibR, capDR, rodR, cpMask, cpRGB}, ...
    'Real Mitochondria — Classical (OGS-smoothed) vs Cellpose', ...
    panelTitles, cmaps);

% =========================================================================
% 8.  Figure 2 — Zoom: cluster 1
% =========================================================================
figure(2);
set(gcf, 'Name', 'Real Image — Zoom cluster 1', 'NumberTitle', 'off', ...
         'Color', 'k', 'Position', [30 340 1820 380]);

mitoFillFigure(2, ...
    {imcrop(Ireal, roi1), imcrop(logR, roi1), imcrop(fibR, roi1), ...
     imcrop(capDR, roi1), imcrop(rodR, roi1), ...
     imcrop(cpMask, roi1), imcrop(cpRGB, roi1)}, ...
    sprintf('Cluster 1  [x=%d  y=%d  %d×%d px]', roi1(1), roi1(2), roi1(3), roi1(4)), ...
    panelTitles, cmaps);

% =========================================================================
% 9.  Figure 3 — Zoom: cluster 2
% =========================================================================
figure(3);
set(gcf, 'Name', 'Real Image — Zoom cluster 2', 'NumberTitle', 'off', ...
         'Color', 'k', 'Position', [30 0 1820 380]);

mitoFillFigure(3, ...
    {imcrop(Ireal, roi2), imcrop(logR, roi2), imcrop(fibR, roi2), ...
     imcrop(capDR, roi2), imcrop(rodR, roi2), ...
     imcrop(cpMask, roi2), imcrop(cpRGB, roi2)}, ...
    sprintf('Cluster 2  [x=%d  y=%d  %d×%d px]', roi2(1), roi2(2), roi2(3), roi2(4)), ...
    panelTitles, cmaps);

% =========================================================================
% 10. Synthetic figures (toggled off)
% =========================================================================
if showSynthetic
    fprintf('\n[showSynthetic = true — no synthetic figures defined in this version]\n');
end

fprintf('\nDone. 3 figures generated.\n');


% =========================================================================
% LOCAL FUNCTION — must be at the end of the script file (R2016b+)
% =========================================================================
function mitoFillFigure(figNum, panels, figTitle, panelTitles, cmaps)
% mitoFillFigure  Populate a 1×N panel figure with black background.
%
%   figNum      - integer figure number (figure must already exist)
%   panels      - 1×N cell array of images
%   figTitle    - sgtitle string
%   panelTitles - 1×N cell of per-panel title strings
%   cmaps       - 1×N cell of colormap names; [] signals truecolor RGB
    figure(figNum);
    N = numel(panels);
    for k = 1:N
        ax = subplot(1, N, k);
        if isempty(cmaps{k})
            imshow(panels{k});          % truecolor RGB (label2rgb output)
        else
            imshow(panels{k}, []);
            colormap(ax, cmaps{k});
        end
        title(panelTitles{k}, 'Color', 'w', 'FontSize', 10);
        set(ax, 'XColor', 'none', 'YColor', 'none');
    end
    sgtitle(figTitle, 'Color', 'w', 'FontSize', 13, 'FontWeight', 'bold');
end

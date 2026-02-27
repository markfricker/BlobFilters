% demoMitoEnhance.m
%
% Demo script: synthetic + real confocal mitochondria images → three enhancers
%
%   logEnhance      – isotropic LoG bank (blob / puncta detection)
%   fiberEnhance    – fibermetric-based tubular enhancer (IPT R2018b+)
%   capsuleEnhance  – oriented capsule bank, single and DoC modes
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

% =========================================================================
% 3.  Run enhancers on BOTH images
% =========================================================================
images  = {Isynth, 'Synthetic'; Ireal, 'Real (mitImage)'};
results = cell(2, 4);   % {logR, fibR, capSR, capDR} per image

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
% 5.  Figure 2 — synthetic image comparison panel
% =========================================================================
titles2 = {'Raw (synthetic)', ...
           'logEnhance', ...
           'fiberEnhance', ...
           'capsule single', ...
           'capsule DoC'};
panels2 = {Isynth, results{1,1}, results{1,2}, results{1,3}, results{1,4}};

figure(2);
set(gcf,'Name','Synthetic — Enhancement Comparison','NumberTitle','off', ...
        'Color','k','Position',[30 30 1250 840]);
cmaps = {'gray','hot','hot','hot','hot'};
for k = 1:5
    ax = subplot(2,3,k);
    imshow(panels2{k},[]); colormap(ax, cmaps{k});
    title(titles2{k},'Color','w','FontSize',10);
    set(ax,'XColor','none','YColor','none');
end
subplot(2,3,6); axis off;
sgtitle('Synthetic Image — Enhancement Comparison', ...
        'Color','w','FontSize',13,'FontWeight','bold');

% =========================================================================
% 6.  Figure 3 — real image comparison panel
% =========================================================================
titles3 = {'Raw (real)', ...
           'logEnhance', ...
           'fiberEnhance', ...
           'capsule single', ...
           'capsule DoC'};
panels3 = {Ireal, results{2,1}, results{2,2}, results{2,3}, results{2,4}};

figure(3);
set(gcf,'Name','Real Image — Enhancement Comparison','NumberTitle','off', ...
        'Color','k','Position',[60 60 1250 840]);
for k = 1:5
    ax = subplot(2,3,k);
    imshow(panels3{k},[]); colormap(ax, cmaps{k});
    title(titles3{k},'Color','w','FontSize',10);
    set(ax,'XColor','none','YColor','none');
end
subplot(2,3,6); axis off;
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
diff_real = results{2,4} - results{2,3};

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

fprintf('\nDone. Five figures generated.\n');
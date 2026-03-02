% exportFigures.m
%
% Generates all figures for the BlobFilters LaTeX manual and saves them as
% PDF files in docs/figures/.  Run this script once before compiling
% BlobFilters_manual.tex.
%
% USAGE
%   cd <repo_root>
%   addpath src demos
%   run docs/exportFigures
%
% OUTPUT FILES  (written to docs/figures/)
%   fig01_synthetic_annotated.pdf   – annotated synthetic test image
%   fig02_comparison_synth.pdf      – 6-panel enhancer comparison (synthetic)
%   fig03_comparison_real.pdf       – 6-panel enhancer comparison (real)
%   fig04_coherence_separation.pdf  – structure tensor coherence (Fig 7)
%   fig05_rod_vs_disk.pdf           – rod vs disk granulometry (Fig 8)
%   fig06_prefilter_raw.pdf         – pre-filter baseline (demoPrefilter Fig 1)
%   fig07_prefilter_pm.pdf          – Perona-Malik (demoPrefilter Fig 2)
%   fig08_prefilter_guided.pdf      – guided filter (demoPrefilter Fig 3)
%   fig09_prefilter_ogs.pdf         – oriented Gauss (demoPrefilter Fig 4)
%   fig10_cellpose.pdf              – Cellpose segmentation (demoMitoEnhance Fig 9)
%                                     OPTIONAL: only written when the Cellpose
%                                     add-on is installed and Fig 9 was created.
%
% NOTE
%   Figures are exported at 300 DPI using exportgraphics (R2020a+).
%   On older releases replace exportgraphics with print -dpdf -r300.

close all;
set(0, 'DefaultFigureVisible', 'off');   % suppress display during export

rootDir  = fullfile(fileparts(mfilename('fullpath')), '..');
outDir   = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(outDir, 'dir'), mkdir(outDir); end

addpath(fullfile(rootDir, 'src'));
addpath(fullfile(rootDir, 'demos'));

% =========================================================================
% 1.  Run demoMitoEnhance (generates Figures 1-8, plus optional Fig 9)
% =========================================================================
fprintf('Running demoMitoEnhance...\n');
demoMitoEnhance;

% demoMitoEnhance calls 'clear' which wipes this script's workspace.
% Recompute outDir and restore headless mode before saving figures.
outDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
set(0, 'DefaultFigureVisible', 'off');

saveFig(figure(1), outDir, 'fig01_synthetic_annotated.pdf');
saveFig(figure(2), outDir, 'fig02_comparison_synth.pdf');
saveFig(figure(3), outDir, 'fig03_comparison_real.pdf');
saveFig(figure(7), outDir, 'fig04_coherence_separation.pdf');
saveFig(figure(8), outDir, 'fig05_rod_vs_disk.pdf');

% Figure 9 (Cellpose) is optional — only created when the add-on is installed.
if ishandle(9)
    saveFig(figure(9), outDir, 'fig10_cellpose.pdf');
end
close all;

% =========================================================================
% 2.  Run demoPrefilter (generates Figures 1-4)
% =========================================================================
fprintf('Running demoPrefilter...\n');
demoPrefilter;

% demoPrefilter also calls 'clear' — recompute outDir again.
outDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
set(0, 'DefaultFigureVisible', 'off');

saveFig(figure(1), outDir, 'fig06_prefilter_raw.pdf');
saveFig(figure(2), outDir, 'fig07_prefilter_pm.pdf');
saveFig(figure(3), outDir, 'fig08_prefilter_guided.pdf');
saveFig(figure(4), outDir, 'fig09_prefilter_ogs.pdf');
close all;

fprintf('\nAll figures exported to %s\n', outDir);

% =========================================================================
% Helper: export a figure handle to PDF, then close it.
%
% outDir is passed explicitly because local functions defined at the bottom
% of a script do not share the script's base workspace in MATLAB.
% =========================================================================
function saveFig(h, outDir, filename)
    outPath = fullfile(outDir, filename);
    exportgraphics(h, outPath, 'ContentType', 'vector', 'Resolution', 300);
    fprintf('  saved %s\n', filename);
    close(h);
end

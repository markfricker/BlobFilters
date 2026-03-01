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
% Helper: export a figure handle to PDF, then close it
% =========================================================================
    function saveFig(h, filename)
        outPath = fullfile(outDir, filename);
        exportgraphics(h, outPath, 'ContentType', 'vector', 'Resolution', 300);
        fprintf('  saved %s\n', filename);
        close(h);
    end

% =========================================================================
% 1.  Run demoMitoEnhance (generates Figures 1-8)
% =========================================================================
fprintf('Running demoMitoEnhance...\n');
demoMitoEnhance;

saveFig(figure(1), 'fig01_synthetic_annotated.pdf');
saveFig(figure(2), 'fig02_comparison_synth.pdf');
saveFig(figure(3), 'fig03_comparison_real.pdf');
saveFig(figure(7), 'fig04_coherence_separation.pdf');
saveFig(figure(8), 'fig05_rod_vs_disk.pdf');
close all;

% =========================================================================
% 2.  Run demoPrefilter (generates Figures 1-4)
% =========================================================================
fprintf('Running demoPrefilter...\n');
demoPrefilter;

saveFig(figure(1), 'fig06_prefilter_raw.pdf');
saveFig(figure(2), 'fig07_prefilter_pm.pdf');
saveFig(figure(3), 'fig08_prefilter_guided.pdf');
saveFig(figure(4), 'fig09_prefilter_ogs.pdf');
close all;

fprintf('\nAll figures exported to %s\n', outDir);

function [R, L] = cellposeEnhance(I, params)
% cellposeEnhance  Deep-learning cell/organelle segmenter via Cellpose
%
% USAGE
%   [R, L] = cellposeEnhance(I)
%   [R, L] = cellposeEnhance(I, params)
%
% INPUTS
%   I       - 2-D grayscale image (any numeric class); converted to
%             single [0,1] internally via im2single.
%   params  - optional struct with fields:
%            .model    = 'cyto3'  % Cellpose model name or absolute path
%                                 % to a custom model file
%            .diameter = 10       % expected object diameter in pixels;
%                                 % set to 0 to trigger auto-estimation
%            .normalize = true    % included for API consistency with other
%                                 % BlobFilters functions; output R is
%                                 % always binary [0,1]
%
% OUTPUTS
%   R       - binary enhancement map, single precision, same size as I.
%             Pixels belonging to detected objects are 1; background is 0.
%             Compatible with the [0,1] convention of all other BlobFilters
%             enhancers and can be used as a segmentation mask directly.
%   L       - label image, uint16, same size as I.
%             Each detected object carries a unique positive integer label.
%             Background = 0.  Use label2rgb(L) for display.
%
% REQUIREMENTS
%   Medical Imaging Toolbox (MathWorks)
%   Medical Imaging Toolbox Interface for Cellpose Library Add-On
%     (install via Home > Add-Ons > search "Cellpose")
%   Cellpose Python package installed in the Python environment that
%     MATLAB's pyenv points to:
%       pip install cellpose
%
% OVERVIEW
%   Wraps the MathWorks cellpose() + segmentCells2D() API introduced in the
%   Medical Imaging Toolbox Interface for Cellpose Library add-on.  A
%   Cellpose model object is created on each call; for batch processing,
%   cache the model object externally and call segmentCells2D() directly to
%   avoid repeated model loading overhead.
%
%   MODEL SELECTION
%   'cyto3'  — Cellpose 3 super-generalist cytoplasm model trained on 9
%              heterogeneous datasets (Nature Methods 2025).  Good default
%              for fluorescence mitochondria images.
%   'cpsam'  — Cellpose-SAM (Cellpose 4); SAM ViT backbone fine-tuned on
%              diverse biological images.  Current best general model.
%   'nuclei' — Nuclear segmentation model.
%   Custom   — Pass an absolute path to a model file trained with the
%              Cellpose GUI (human-in-the-loop mode) for best mitochondria
%              specificity after fine-tuning on ~100 annotated images.
%
%   DIAMETER PARAMETER
%   Set diameter to the expected object cross-section in pixels (not the
%   long axis).  For ~8 px wide mitochondria, diameter = 8-12 is
%   appropriate.  Setting diameter = 0 triggers Cellpose's built-in
%   automatic size estimation via a companion size model (slower, ~2x).
%
%   OUTPUT CONVENTION
%   Unlike all other BlobFilters enhancers, which return continuous
%   probability-like maps, cellposeEnhance returns a binary map (R) and
%   an integer label image (L).  The binary map R can be used in the same
%   pipeline positions as other enhancer outputs (thresholding, logical
%   masking), while L additionally encodes individual object identity.
%
% NOTES
%   - The function checks for the cellpose() function at runtime and raises
%     an informative error if the add-on is not installed.
%   - Cellpose internally normalises the input image; passing normalised
%     single [0,1] or raw uint8/uint16 both give equivalent results.
%   - For 3-D stacks use segmentCells3D() from the Medical Imaging Toolbox
%     directly; cellposeEnhance processes 2-D slices only.
%
% REFERENCES
%   Stringer~C., Wang~T., Michaelos~M. \& Pachitariu~M. (2021)
%   Cellpose: a generalist algorithm for cellular segmentation.
%   Nature Methods 18:100-106.
%   https://doi.org/10.1038/s41592-020-01018-x
%     -> Original Cellpose paper and cyto/nuclei models.
%
%   Pachitariu~M. \& Stringer~C. (2022)
%   Cellpose 2.0: how to train your own model.
%   Nature Methods 19:1634-1641.
%   https://doi.org/10.1038/s41592-022-01663-4
%     -> Human-in-the-loop fine-tuning for custom models.
%
%   Stringer~C. \& Pachitariu~M. (2025)
%   Cellpose3: one-click image restoration for improved cellular
%   segmentation. Nature Methods 22:592-599.
%   https://doi.org/10.1038/s41592-025-02595-5
%     -> cyto3 model and Cellpose3 paper.
%
%   MathWorks (2025) Medical Imaging Toolbox Interface for Cellpose Library.
%   https://www.mathworks.com/help/medical-imaging/ref/cellpose.html
%     -> MATLAB wrapper used by this function.
%
% EXAMPLE
%   % Basic segmentation with default cyto3 model
%   pCP.model    = 'cyto3';
%   pCP.diameter = 10;       % ~8 px wide mitochondria
%   [R, L] = cellposeEnhance(I, pCP);
%
%   % Display results
%   figure;
%   subplot(1,3,1); imshow(I,[]);          title('Raw');
%   subplot(1,3,2); imshow(R,[]);          title('Binary mask');
%   subplot(1,3,3); imshow(label2rgb(L,'hsv','k')); title('Labels');
%
%   % Use binary mask as weight on another enhancer
%   R_rod = rodGranulometryEnhance(I) .* R;
%
%   % Custom fine-tuned model
%   pCP.model = 'C:\models\mito_finetuned';
%   [R, L] = cellposeEnhance(I, pCP);
%
% See also: fiberEnhance, rodGranulometryEnhance, structureTensorEnhance

% --- defaults ---------------------------------------------------------------
if nargin < 2, params = struct(); end
if ~isfield(params, 'model'),     params.model     = 'cyto3'; end
if ~isfield(params, 'diameter'),  params.diameter  = 10;      end
if ~isfield(params, 'normalize'), params.normalize = true;    end

% --- input validation -------------------------------------------------------
if size(I, 3) > 1
    error('cellposeEnhance:badInput', ...
          'cellposeEnhance: expected 2-D grayscale image, got %d-channel input.', ...
          size(I,3));
end
if ~exist('cellpose', 'file')
    error('cellposeEnhance:notFound', ...
          ['cellposeEnhance: cellpose() not found.\n' ...
           'Install the "Medical Imaging Toolbox Interface for Cellpose Library"\n' ...
           'via Home > Add-Ons, then ensure Cellpose is installed in your\n' ...
           'Python environment (pip install cellpose).']);
end

I = im2single(I);

% --- run Cellpose -----------------------------------------------------------
cp = cellpose(Model=params.model);
L  = uint16(segmentCells2D(cp, I, ImageCellDiameter=params.diameter));

% Binary enhancement map: 1 inside detected objects, 0 in background.
% Consistent with the [0,1] output convention of all other BlobFilters
% enhancers.
R = single(L > 0);
end

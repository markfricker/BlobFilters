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
%            .model         = 'cyto3' % Cellpose model name or absolute path
%                                     % to a custom model file
%            .diameter      = 10      % expected object diameter in pixels
%            .cellProb      = 0       % cell probability threshold in [−6,6];
%                                     % lower values detect fainter / smaller
%                                     % objects.  0 found optimal for
%                                     % mitochondria fluorescence (sweepCellpose
%                                     % 2026-03-03); try −1 for slightly higher
%                                     % recall at cost of more false positives.
%            .flowThreshold = 0.4     % flow-error threshold in [0.1,3];
%                                     % higher values detect more objects
%                                     % at the cost of boundary precision
%            .nIter         = 0       % Cellpose gradient-descent iterations.
%                                     % 0 = use Cellpose default (~200).
%                                     % Set to 2000 for cyto3 on long thin
%                                     % structures (mitochondria, bacteria).
%                                     % Uses the Python Cellpose API directly
%                                     % (MathWorks segmentCells2D does not
%                                     % expose this parameter).
%            .normalize     = true    % API consistency; R is always binary
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
%   'cyto3'  — Cellpose 3 super-generalist cytoplasm model (DEFAULT).
%              Trained on 9 heterogeneous datasets (Nature Methods 2025);
%              best general-purpose choice for fluorescence mito images.
%              Requires Python cellpose >= 3.x (the MATLAB add-on's
%              downloadCellposeModels does not include cyto3; the function
%              downloads it automatically via Python on first call).
%   'bact_fluor_cp3' — Cellpose 3 model trained specifically on fluorescent
%              bacterial images.  Bacteria and mitochondria share similar
%              elongated morphology, making this model a strong alternative
%              to cyto3 for organelle segmentation.  Weights are downloaded
%              automatically on first call (~80 MB).  Use nIter = 2000 for
%              best results on long thin mitochondria.
%   'cyto2'  — Cellpose 2 cytoplasm model.  Weights are in the MATLAB
%              add-on's registry (downloadCellposeModels); loaded by file
%              path, no Python download needed.  Safe fallback when
%              Python cellpose 4.x is installed (cyto3 → cpsam in 4.x).
%   'cpsam'  — Cellpose-SAM (Cellpose 4 only); SAM ViT backbone fine-tuned
%              on diverse biological cell images.  Not suitable for small
%              organelles such as mitochondria (returns blank masks with
%              ImageCellDiameter = 10–80).
%   'nuclei' — Nuclear segmentation model.
%   Custom   — Pass an absolute path to a model file trained with the
%              Cellpose GUI (human-in-the-loop mode) for best mitochondria
%              specificity after fine-tuning on ~100 annotated images.
%
%   DIAMETER PARAMETER
%   Set diameter to the expected object cross-section in pixels (not the
%   long axis).  For ~8 px wide mitochondria, diameter = 8-12 is
%   appropriate.  Cellpose internally rescales the image so that objects
%   appear ~30 px wide before inference, so diameter accuracy matters more
%   for boundary quality than for whether objects are detected at all.
%
%   SENSITIVITY PARAMETERS
%   The two main knobs for improving recall on faint or small objects:
%
%   cellProb (CellThreshold, default 0, range −6 to 6):
%     Lower values accept lower-probability detections.  Empirically,
%     0 gives the best boundary quality for mitochondria fluorescence
%     images (sweepCellpose, 2026-03-03); try −1 for slightly higher
%     recall.  Values below −3 tend to produce excessive false positives.
%
%   flowThreshold (FlowErrorThreshold, default 0.4, range 0.1 to 3):
%     Acceptable error between predicted and reconstructed flow fields.
%     Increase to 0.8–1.0 to recover fragmentary or elongated objects
%     whose flow fields are imperfect; at the cost of rougher boundaries.
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
%   - When nIter > 0, cellposeEnhance bypasses segmentCells2D entirely and
%     calls the Python CellposeModel.eval() API directly (since the MATLAB
%     wrapper does not expose niter).  This requires MATLAB R2022a+ for
%     correct numpy <-> MATLAB array conversion.  In this path the model is
%     resolved by Python directly, so cpResolveModelPath is not called.
%   - The MATLAB add-on's model registry (downloadCellposeModels) only
%     covers models up to Cellpose 2 ('cyto', 'cyto2', 'nuclei', ...).
%     Newer models ('cyto3', 'bact_fluor_cp3', 'cpsam') are NOT in that
%     list and must be referenced by absolute file path.
%     cellposeEnhance resolves this automatically: if cellpose() rejects a
%     model name, cpResolveModelPath uses the Python cellpose package to
%     download the weights (first call only) and returns the cached file
%     path.  When nIter > 0 the Python path handles model resolution
%     natively without cpResolveModelPath.
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
%   % Increase recall for faint/small mitochondria
%   pCP.cellProb      = -2;  % accept lower-probability detections
%   pCP.flowThreshold = 0.8; % tolerate imperfect flow fields
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
%   % cyto3 with increased niter for elongated mitochondria (literature rec.)
%   pCP.model = 'cyto3';
%   pCP.nIter = 2000;    % uses Python API directly; not exposed by MATLAB add-on
%   [R, L] = cellposeEnhance(I, pCP);
%
%   % bact_fluor_cp3: Cellpose3 model trained on fluorescent bacteria (good for mito)
%   pCP.model = 'bact_fluor_cp3';
%   pCP.nIter = 2000;
%   [R, L] = cellposeEnhance(I, pCP);
%
%   % Custom fine-tuned model
%   pCP.model = 'C:\models\mito_finetuned';
%   [R, L] = cellposeEnhance(I, pCP);
%
% See also: fiberEnhance, rodGranulometryEnhance, structureTensorEnhance

% --- defaults ---------------------------------------------------------------
if nargin < 2, params = struct(); end
if ~isfield(params, 'model'),         params.model         = 'cyto3'; end
if ~isfield(params, 'diameter'),      params.diameter      = 10;      end
if ~isfield(params, 'cellProb'),      params.cellProb      = -2;       end
if ~isfield(params, 'flowThreshold'), params.flowThreshold = 0.8;     end
if ~isfield(params, 'nIter'),         params.nIter         = 0;       end
if ~isfield(params, 'normalize'),     params.normalize     = true;    end

% --- input validation -------------------------------------------------------
if size(I, 3) > 1
    error('cellposeEnhance:badInput', ...
          'cellposeEnhance: expected 2-D grayscale image, got %d-channel input.', ...
          size(I,3));
end
if exist('cellpose', 'file') == 0
    error('cellposeEnhance:notFound', ...
          ['cellposeEnhance: cellpose() not found.\n' ...
           'Install the "Medical Imaging Toolbox Interface for Cellpose Library"\n' ...
           'via Home > Add-Ons, then ensure Cellpose is installed in your\n' ...
           'Python environment (pip install cellpose).']);
end

I = im2single(I);

% --- nIter > 0: bypass segmentCells2D and call Python directly --------------
% MathWorks segmentCells2D does not expose the niter parameter.
% When the user sets nIter, we call CellposeModel.eval() via pyrun instead.
if params.nIter > 0
    [R, L] = cpRunWithNIter(I, params.model, params.diameter, ...
                            params.cellProb, params.flowThreshold, ...
                            params.nIter);
    return
end

% --- run Cellpose -----------------------------------------------------------
% The MATLAB add-on only recognises models in its own registry (up to
% Cellpose 2).  Newer names such as 'cyto3' / 'cpsam' require an absolute
% file path.  If the named lookup fails, fall back to Python resolution.
try
    cp = cellpose(Model=params.model);
catch ME
    if contains(ME.message, 'Unable to find model file')
        params.model = cpResolveModelPath(params.model);
        cp = cellpose(Model=params.model);
    else
        rethrow(ME);
    end
end
L  = uint16(segmentCells2D(cp, I, ...
         ImageCellDiameter = params.diameter, ...
         CellThreshold     = params.cellProb, ...
         FlowErrorThreshold= params.flowThreshold));

% Binary enhancement map: 1 inside detected objects, 0 in background.
% Consistent with the [0,1] output convention of all other BlobFilters
% enhancers.
R = single(L > 0);
end

% ---- local helper: Python Cellpose with niter support ----------------------
function [R, L] = cpRunWithNIter(I, modelName, diameter, cellProb, flowThreshold, nIter)
% cpRunWithNIter  Cellpose segmentation via direct Python call.
%   Supports the niter parameter, which MathWorks segmentCells2D() does not
%   expose.  Named models ('cyto3', 'bact_fluor_cp3', etc.) are resolved by
%   the Python Cellpose package directly (downloading weights on first call).
%   Custom paths are passed via pretrained_model.  Requires MATLAB R2022a+.
%
%   I is already im2single [0,1] on entry.

% Select the Python constructor keyword based on whether modelName is a path
if contains(modelName, filesep) || contains(modelName, '/')
    ctorStr = "CellposeModel(pretrained_model='" + modelName + "')";
else
    ctorStr = "CellposeModel(model_type='" + modelName + "')";
end

masks = pyrun([ ...
    "from cellpose.models import CellposeModel; " ...
    "import numpy as np; " ...
    "m = " + ctorStr + "; " ...
    "msk, _, _ = m.eval([img.astype(np.float32)], " ...
    "    diameter=float(diam), niter=int(nit), " ...
    "    cellprob_threshold=float(cp), " ...
    "    flow_threshold=float(ft), do_3D=False); " ...
    "out = msk[0].astype(np.uint16)"], ...
    "out", img=double(I), diam=diameter, nit=int32(nIter), ...
    cp=cellProb, ft=flowThreshold);

% Convert numpy uint16 [H x W] → MATLAB uint16 matrix.
% double() intermediate handles numpy→MATLAB conversion (R2022a+).
L = uint16(double(masks));
if ~isequal(size(L), size(I))
    L = L';   % transpose if MATLAB/Python axis order is swapped
end
R = single(L > 0);
end


% ---- local helper -----------------------------------------------------------
function p = cpResolveModelPath(modelName)
% cpResolveModelPath  Resolve a Cellpose model name to its cached file path.
%
%   Uses the Python cellpose package to instantiate the requested model
%   (downloading weights if they are not already cached) and returns the
%   absolute path to the weight file.  This path can be passed directly to
%   the MATLAB cellpose() constructor as a custom model.
%
%   Handles the Cellpose 4+ behaviour where 'cyto3' is silently redirected
%   to 'cpsam' — the returned path reflects whichever file was actually
%   cached, not the name requested.
p = char(pyrun( ...
    "from cellpose.models import CellposeModel; " + ...
    "m = CellposeModel(model_type='" + modelName + "'); " + ...
    "out = m.pretrained_model[0] if isinstance(m.pretrained_model, list) " + ...
    "else m.pretrained_model", ...
    "out"));
end

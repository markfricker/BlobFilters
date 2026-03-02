classdef testBlobFilters < matlab.unittest.TestCase
% testBlobFilters  Unit tests for the BlobFilters toolbox
%
% USAGE
%   Run all tests from the MATLAB command window:
%       results = runtests('tests/testBlobFilters');
%       table(results)
%
%   Run a single test method:
%       runtests('tests/testBlobFilters/testLog_blobResponse')
%
% COVERAGE
%   logEnhance            – 7 tests
%   fiberEnhance          – 7 tests
%   capsuleEnhance        – 7 tests
%   granulometryEnhance   – 7 tests
%   rodGranulometryEnhance– 7 tests
%   structureTensorEnhance– 8 tests
%   orientedGaussSmooth   – 6 tests
%
% REQUIREMENTS
%   MATLAB R2013b+ (matlab.unittest framework)
%   Image Processing Toolbox (imfilter, imgaussfilt, fibermetric R2018b+)
%   All src/ functions on the path (added automatically by TestClassSetup).

    properties (Constant)
        Tol = single(1e-4);   % absolute tolerance for floating-point checks
    end

    % =====================================================================
    % Path setup – runs once before any test in the class
    % =====================================================================
    methods (TestClassSetup)
        function addSrcPath(tc)
            rootDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(rootDir, 'src'));
            addpath(fullfile(rootDir, 'demos'));
        end
    end

    % =====================================================================
    % Shared synthetic test images
    %
    %   uniformImage  – constant grey (no structure)
    %   blobImage     – centred isotropic Gaussian (blob / puncta)
    %   rodImage      – short horizontal bright rod (16 px, Gaussian x-section)
    %                   Rod is shorter than the default max length in
    %                   rodGranulometryEnhance so the pattern spectrum has a
    %                   well-defined peak at the rod scale.
    % =====================================================================
    methods (Static, Access = private)

        function I = uniformImage()
            I = 0.5 * ones(64, 64, 'single');
        end

        function I = blobImage()
            [xx, yy] = meshgrid(1:64, 1:64);
            I = single(exp(-((xx-32).^2 + (yy-32).^2) / (2*5^2)));
        end

        function I = rodImage()
            % Horizontal rod: Gaussian cross-section (sigma_y=3 px), 16 px long.
            % Centred at (row=32, col=32).  Cols 24-39 are lit.
            [xx, yy] = meshgrid(1:64, 1:64);
            I = single(exp(-(yy - 32).^2 / (2*3^2)) .* (xx >= 24 & xx <= 39));
            I = I / max(I(:));
        end

    end

    % =====================================================================
    %  logEnhance
    % =====================================================================
    methods (Test)

        function testLog_smoke(tc)
            R = logEnhance(testBlobFilters.blobImage());
            tc.verifyNotEmpty(R);
        end

        function testLog_outputSize(tc)
            I = testBlobFilters.blobImage();
            tc.verifySize(logEnhance(I), size(I));
        end

        function testLog_outputClass(tc)
            tc.verifyClass(logEnhance(testBlobFilters.blobImage()), 'single');
        end

        function testLog_normalizedRange(tc)
            R = logEnhance(testBlobFilters.blobImage(), struct('normalize', true));
            tc.verifyGreaterThanOrEqual(min(R(:)), 0);
            tc.verifyLessThanOrEqual(max(R(:)), 1 + tc.Tol);
        end

        function testLog_normalizedMax(tc)
            R = logEnhance(testBlobFilters.blobImage(), struct('normalize', true));
            tc.verifyEqual(max(R(:)), single(1), 'AbsTol', tc.Tol);
        end

        function testLog_uniformImage_givesZero(tc)
            % Uniform image has zero Laplacian -> zero LoG response
            R = logEnhance(testBlobFilters.uniformImage(), struct('normalize', false));
            tc.verifyLessThanOrEqual(max(abs(R(:))), tc.Tol);
        end

        function testLog_blobResponse(tc)
            % A centred Gaussian blob should give maximum response at centre
            I = testBlobFilters.blobImage();
            R = logEnhance(I, struct('sigmas', [3 4 5 6], 'normalize', true));
            tc.verifyGreaterThan(R(32, 32), 0.5);
        end

    end

    % =====================================================================
    %  fiberEnhance
    % =====================================================================
    methods (Test)

        function testFiber_smoke(tc)
            tc.verifyNotEmpty(fiberEnhance(testBlobFilters.rodImage()));
        end

        function testFiber_outputSize(tc)
            I = testBlobFilters.rodImage();
            tc.verifySize(fiberEnhance(I), size(I));
        end

        function testFiber_outputClass(tc)
            tc.verifyClass(fiberEnhance(testBlobFilters.rodImage()), 'single');
        end

        function testFiber_normalizedRange(tc)
            R = fiberEnhance(testBlobFilters.rodImage(), struct('normalize', true));
            tc.verifyGreaterThanOrEqual(min(R(:)), 0);
            tc.verifyLessThanOrEqual(max(R(:)), 1 + tc.Tol);
        end

        function testFiber_uniformImage_givesZero(tc)
            R = fiberEnhance(testBlobFilters.uniformImage(), struct('normalize', false));
            tc.verifyLessThanOrEqual(max(abs(R(:))), tc.Tol);
        end

        function testFiber_stackMode(tc)
            p.widths = [4 6]; p.multimode = 'stack'; p.normalize = true;
            R = fiberEnhance(testBlobFilters.rodImage(), p);
            tc.verifySize(R, size(testBlobFilters.rodImage()));
        end

        function testFiber_rodResponse(tc)
            % Hessian-based tubularity should respond strongly on a bright rod
            p.widths = [4 5 6]; p.multimode = 'builtin'; p.normalize = true;
            R = fiberEnhance(testBlobFilters.rodImage(), p);
            rodMean = mean(R(32, 24:39));
            bgMean  = mean(mean(R(1:10, 1:10)));
            tc.verifyGreaterThan(rodMean, bgMean);
        end

    end

    % =====================================================================
    %  capsuleEnhance
    % =====================================================================
    methods (Test)

        function testCapsule_smoke(tc)
            tc.verifyNotEmpty(capsuleEnhance(testBlobFilters.rodImage()));
        end

        function testCapsule_outputSize(tc)
            I = testBlobFilters.rodImage();
            tc.verifySize(capsuleEnhance(I), size(I));
        end

        function testCapsule_outputClass(tc)
            tc.verifyClass(capsuleEnhance(testBlobFilters.rodImage()), 'single');
        end

        function testCapsule_normalizedRange(tc)
            R = capsuleEnhance(testBlobFilters.rodImage(), struct('normalize', true));
            tc.verifyGreaterThanOrEqual(min(R(:)), 0);
            tc.verifyLessThanOrEqual(max(R(:)), 1 + tc.Tol);
        end

        function testCapsule_singleMode(tc)
            p.mode = 'single'; p.lengths = [12 16]; p.width = 6;
            p.orientations = 4; p.normalize = true;
            R = capsuleEnhance(testBlobFilters.rodImage(), p);
            tc.verifySize(R, size(testBlobFilters.rodImage()));
        end

        function testCapsule_docMode(tc)
            p.mode = 'doc'; p.lengths = [12 16]; p.width = 6;
            p.wideWidth = 12; p.orientations = 4; p.normalize = true;
            R = capsuleEnhance(testBlobFilters.rodImage(), p);
            tc.verifySize(R, size(testBlobFilters.rodImage()));
        end

        function testCapsule_rodResponse(tc)
            % Capsule at 0 deg should respond along the horizontal rod
            p.mode = 'single'; p.lengths = [12 16]; p.width = 6;
            p.orientations = 8; p.normalize = true;
            R = capsuleEnhance(testBlobFilters.rodImage(), p);
            rodMean = mean(R(32, 24:39));
            bgMean  = mean(mean(R(1:10, 1:10)));
            tc.verifyGreaterThan(rodMean, bgMean);
        end

    end

    % =====================================================================
    %  granulometryEnhance
    % =====================================================================
    methods (Test)

        function testGran_smoke(tc)
            tc.verifyNotEmpty(granulometryEnhance(testBlobFilters.blobImage()));
        end

        function testGran_outputSize(tc)
            I = testBlobFilters.blobImage();
            tc.verifySize(granulometryEnhance(I), size(I));
        end

        function testGran_outputClass(tc)
            tc.verifyClass(granulometryEnhance(testBlobFilters.blobImage()), 'single');
        end

        function testGran_normalizedRange(tc)
            G = granulometryEnhance(testBlobFilters.blobImage(), struct('normalize', true));
            tc.verifyGreaterThanOrEqual(min(G(:)), 0);
            tc.verifyLessThanOrEqual(max(G(:)), 1 + tc.Tol);
        end

        function testGran_uniformImage_givesZero(tc)
            G = granulometryEnhance(testBlobFilters.uniformImage(), struct('normalize', false));
            tc.verifyLessThanOrEqual(max(abs(G(:))), tc.Tol);
        end

        function testGran_tooFewScales_warns(tc)
            tc.verifyWarning( ...
                @() granulometryEnhance(testBlobFilters.blobImage(), struct('sigmas', 4)), ...
                'granulometryEnhance:tooFewScales');
        end

        function testGran_blobResponse(tc)
            % Disk opening should respond at the blob scale
            p.sigmas = [2 4 6 8]; p.normalize = true;
            G = granulometryEnhance(testBlobFilters.blobImage(), p);
            blobMean = mean(G(30:34, 30:34), 'all');
            bgMean   = mean(G(1:8,   1:8  ), 'all');
            tc.verifyGreaterThan(blobMean, bgMean);
        end

    end

    % =====================================================================
    %  rodGranulometryEnhance
    % =====================================================================
    methods (Test)

        function testRodGran_smoke(tc)
            tc.verifyNotEmpty(rodGranulometryEnhance(testBlobFilters.rodImage()));
        end

        function testRodGran_outputSize(tc)
            I = testBlobFilters.rodImage();
            tc.verifySize(rodGranulometryEnhance(I), size(I));
        end

        function testRodGran_outputClass(tc)
            tc.verifyClass(rodGranulometryEnhance(testBlobFilters.rodImage()), 'single');
        end

        function testRodGran_normalizedRange(tc)
            R = rodGranulometryEnhance(testBlobFilters.rodImage(), struct('normalize', true));
            tc.verifyGreaterThanOrEqual(min(R(:)), 0);
            tc.verifyLessThanOrEqual(max(R(:)), 1 + tc.Tol);
        end

        function testRodGran_uniformImage_givesZero(tc)
            R = rodGranulometryEnhance(testBlobFilters.uniformImage(), struct('normalize', false));
            tc.verifyLessThanOrEqual(max(abs(R(:))), tc.Tol);
        end

        function testRodGran_tooFewLengths_warns(tc)
            p.lengths = [16]; % single length -> no residues
            tc.verifyWarning( ...
                @() rodGranulometryEnhance(testBlobFilters.rodImage(), p), ...
                'rodGranulometryEnhance:tooFewLengths');
        end

        function testRodGran_rodSelectivity(tc)
            % Rod (16 px) is shorter than max length in bank (24 px), so
            % pattern spectrum residue O_max(12)-O_max(24) is large at rod.
            p.lengths = [8 12 24]; p.orientations = 8; p.normalize = true;
            R = rodGranulometryEnhance(testBlobFilters.rodImage(), p);
            rodMean = mean(R(32, 24:39));
            bgMean  = mean(mean(R(1:10, 1:10)));
            tc.verifyGreaterThan(rodMean, bgMean);
        end

    end

    % =====================================================================
    %  structureTensorEnhance
    % =====================================================================
    methods (Test)

        function testST_smoke(tc)
            tc.verifyNotEmpty(structureTensorEnhance(testBlobFilters.rodImage()));
        end

        function testST_outputSize(tc)
            I = testBlobFilters.rodImage();
            tc.verifySize(structureTensorEnhance(I), size(I));
        end

        function testST_outputClass(tc)
            tc.verifyClass(structureTensorEnhance(testBlobFilters.rodImage()), 'single');
        end

        function testST_coherenceRange(tc)
            C = structureTensorEnhance(testBlobFilters.rodImage());
            tc.verifyGreaterThanOrEqual(min(C(:)), 0);
            tc.verifyLessThanOrEqual(max(C(:)), 1 + tc.Tol);
        end

        function testST_thetaOutputSize(tc)
            I = testBlobFilters.rodImage();
            [C, theta] = structureTensorEnhance(I);
            tc.verifySize(theta, size(C));
        end

        function testST_thetaRange(tc)
            [~, theta] = structureTensorEnhance(testBlobFilters.rodImage());
            tc.verifyGreaterThanOrEqual(min(theta(:)), single(-pi/2) - tc.Tol);
            tc.verifyLessThanOrEqual(   max(theta(:)), single( pi/2) + tc.Tol);
        end

        function testST_uniformCoherence(tc)
            % No gradients on a uniform image -> coherence should be zero
            C = structureTensorEnhance(testBlobFilters.uniformImage());
            tc.verifyLessThanOrEqual(max(C(:)), tc.Tol);
        end

        function testST_rodOrientation(tc)
            % Horizontal rod -> along-fibre theta near 0 at centre
            p.sigmaGrad = 1; p.sigmaInt = 3; p.normalize = true;
            [~, theta] = structureTensorEnhance(testBlobFilters.rodImage(), p);
            tc.verifyLessThan(abs(theta(32, 32)), single(0.3));  % within ~17 deg
        end

    end

    % =====================================================================
    %  orientedGaussSmooth
    % =====================================================================
    methods (Test)

        function testOGS_smoke(tc)
            tc.verifyNotEmpty(orientedGaussSmooth(testBlobFilters.rodImage()));
        end

        function testOGS_outputSize(tc)
            I = testBlobFilters.rodImage();
            tc.verifySize(orientedGaussSmooth(I), size(I));
        end

        function testOGS_outputClass(tc)
            tc.verifyClass(orientedGaussSmooth(testBlobFilters.rodImage()), 'single');
        end

        function testOGS_outputWithinInputBounds(tc)
            % Weighted average cannot exceed input range
            I   = testBlobFilters.rodImage();
            Ism = orientedGaussSmooth(I);
            tc.verifyGreaterThanOrEqual(min(Ism(:)), min(I(:)) - tc.Tol);
            tc.verifyLessThanOrEqual(   max(Ism(:)), max(I(:)) + tc.Tol);
        end

        function testOGS_uniformPreservation(tc)
            % All oriented Gaussians applied to a constant return that constant;
            % the weighted blend must also return the same constant.
            I   = testBlobFilters.uniformImage();
            Ism = orientedGaussSmooth(I);
            tc.verifyEqual(Ism, I, 'AbsTol', tc.Tol);
        end

        function testOGS_sigmaOrderWarning(tc)
            p.sigmaAlong = 1; p.sigmaAcross = 4;  % inverted: across > along
            tc.verifyWarning( ...
                @() orientedGaussSmooth(testBlobFilters.rodImage(), p), ...
                'orientedGaussSmooth:sigmaOrder');
        end

    end

    % =====================================================================
    %  cellposeEnhance
    %
    %  Tests that require Cellpose use tc.assumeTrue to mark themselves as
    %  Incomplete (not Failed) when the add-on is absent, consistent with
    %  how fiberEnhance tests handle a missing Image Processing Toolbox.
    % =====================================================================
    methods (Test)

        function testCP_smoke(tc)
            tc.assumeTrue(exist('cellpose','file') ~= 0, ...
                'Cellpose add-on not installed — test skipped');
            [R, L] = cellposeEnhance(testBlobFilters.blobImage());
            tc.verifyNotEmpty(R);
            tc.verifyNotEmpty(L);
        end

        function testCP_outputSize(tc)
            tc.assumeTrue(exist('cellpose','file') ~= 0, ...
                'Cellpose add-on not installed — test skipped');
            I = testBlobFilters.blobImage();
            [R, L] = cellposeEnhance(I);
            tc.verifySize(R, size(I));
            tc.verifySize(L, size(I));
        end

        function testCP_outputClass(tc)
            tc.assumeTrue(exist('cellpose','file') ~= 0, ...
                'Cellpose add-on not installed — test skipped');
            [R, L] = cellposeEnhance(testBlobFilters.blobImage());
            tc.verifyClass(R, 'single');
            tc.verifyClass(L, 'uint16');
        end

        function testCP_binaryRange(tc)
            % R must be binary: all values exactly 0 or 1
            tc.assumeTrue(exist('cellpose','file') ~= 0, ...
                'Cellpose add-on not installed — test skipped');
            R = cellposeEnhance(testBlobFilters.blobImage());
            tc.verifyGreaterThanOrEqual(min(R(:)), single(0));
            tc.verifyLessThanOrEqual(   max(R(:)), single(1));
            uniqueVals = unique(R(:));
            tc.verifyTrue(all(uniqueVals == 0 | uniqueVals == 1));
        end

        function testCP_uniformImage_givesZero(tc)
            % A flat uniform image has no cell-like structure; Cellpose
            % should detect no objects and return an all-zero binary mask.
            tc.assumeTrue(exist('cellpose','file') ~= 0, ...
                'Cellpose add-on not installed — test skipped');
            R = cellposeEnhance(testBlobFilters.uniformImage());
            tc.verifyEqual(max(R(:)), single(0));
        end

        function testCP_toolboxMissing_errors(tc)
            % Verify that the function raises a specific, informative error
            % when the Cellpose add-on is not installed.
            % This test is only meaningful when the add-on IS absent;
            % skip it (Incomplete) when Cellpose is installed.
            tc.assumeFalse(exist('cellpose','file') ~= 0, ...
                'Cellpose IS installed — missing-toolbox test not applicable');
            tc.verifyError( ...
                @() cellposeEnhance(testBlobFilters.blobImage()), ...
                'cellposeEnhance:notFound');
        end

    end

end

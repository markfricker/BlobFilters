function chromosomes = segmentMetaphaseChromosomes(img)

%% --- 0. Ensure grayscale ---
if size(img, 3) == 3
    img = rgb2gray(img);
end
img = im2double(img);

%% --- 1. Denoise FIRST (before any enhancement) ---
% Non-local means: excellent at preserving sharp chromosome edges
% while killing the flat noise texture you see in the background.
% Tune DegreeOfSmoothing: 0.005–0.02 (higher = more smoothing)
img_denoised = imnlmfilt(img, 'DegreeOfSmoothing', 0.01);

%% --- 2. Background Correction on denoised image ---
se_bg = strel('disk', 25);
bg = imopen(img_denoised, se_bg);
img_corr = img_denoised - bg;
img_corr = mat2gray(img_corr);

%% --- 3. CLAHE AFTER denoising (now it enhances signal, not noise) ---
img_enh = adapthisteq(img_corr, ...
    'NumTiles',    [6 6], ...
    'ClipLimit',   0.03, ...
    'Distribution','rayleigh');

%% --- 4. DoG blob enhancement (tuned to chromosome size) ---
% Sigma1 should be ~half the chromosome width: 8px/2 = 4
% Sigma2 should be ~chromosome width: 8px
g1 = imgaussfilt(img_enh, 2.0);
g2 = imgaussfilt(img_enh, 6.0);
img_dog = g1 - 0.7 * g2;
img_dog = max(img_dog, 0);        % suppress negative values (background)
img_dog = mat2gray(img_dog);

%% --- 5. Threshold: Otsu on DoG, with conservative multiplier ---
level = graythresh(img_dog);

% Plot histogram so you can tune the multiplier visually
figure('Name','Threshold Diagnostic');
imhist(img_dog); hold on;
xline(level,      'r--', 'LineWidth', 2, 'Label', 'Otsu');
xline(level*0.80, 'g--', 'LineWidth', 2, 'Label', 'Used (×0.80)');
title('DoG Histogram — tune multiplier until threshold sits in valley');

bw = imbinarize(img_dog, level * 0.80);  % raise from 0.65 → 0.80

%% --- 6. Morphological Cleanup ---
bw = bwareaopen(bw, 30);           % remove small specks

se_close = strel('disk', 2);       % smaller closing (less merging)
bw = imclose(bw, se_close);
bw = imfill(bw, 'holes');
bw = bwareaopen(bw, 40);           % second pass after closing

%% --- 7. Watershed separation ---
D       = bwdist(~bw);
D_s     = imgaussfilt(D, 1.5);
D_neg   = imhmin(-D_s, 2.0);       % h=2.0: tune up if over-splitting
L       = watershed(D_neg);
bw_ws   = bw;
bw_ws(L == 0) = 0;

%% --- 8. Connected component filtering ---
CC    = bwconncomp(bw_ws);
props = regionprops(CC, img, ...
    'Area', 'BoundingBox', 'MajorAxisLength', 'MinorAxisLength', ...
    'Solidity', 'Image', 'Extent');

% Chromosome dimensions: width~8px, length 12-40px
% Area range: 8*12=96 to 8*40=320, allow 2x margin each side
minArea     = 60;
maxArea     = 700;
minLength   = 10;
maxLength   = 90;
minSolidity = 0.50;
maxWidth    = 25;    % rejects filamentary noise

area   = [props.Area];
majLen = [props.MajorAxisLength];
minWid = [props.MinorAxisLength];
solid  = [props.Solidity];

% keep = area   >= minArea    & ...
%        area   <= maxArea    & ...
%        majLen >= minLength  & ...
%        majLen <= maxLength  & ...
%        minWid <= maxWidth   & ...
%        solid  >= minSolidity;
keep = true(1,CC.NumObjects);

props_kept = props(keep);
kept_idx   = find(keep);

CC_kept = CC;
CC_kept.NumObjects   = sum(keep);
CC_kept.PixelIdxList = CC.PixelIdxList(kept_idx);

%% --- 9. Package output ---
chromosomes = struct([]);
for k = 1:numel(props_kept)
    bb = props_kept(k).BoundingBox;
    x1 = max(round(bb(1)), 1);
    y1 = max(round(bb(2)), 1);
    x2 = min(round(bb(1)+bb(3)-1), size(img,2));
    y2 = min(round(bb(2)+bb(4)-1), size(img,1));

    chromosomes(k).label    = k;
    chromosomes(k).bbox     = [x1 y1 x2 y2];
    chromosomes(k).mask     = props_kept(k).Image;
    chromosomes(k).image    = img_corr(y1:y2, x1:x2);
    chromosomes(k).area     = props_kept(k).Area;
    chromosomes(k).length   = props_kept(k).MajorAxisLength;
    chromosomes(k).width    = props_kept(k).MinorAxisLength;
    chromosomes(k).solidity = props_kept(k).Solidity;
end

%% --- 10. Visualise ---
figure('Position', [100 100 1400 420]);

subplot(1,4,1);
imshow(img,[]); title('Original');

subplot(1,4,2);
imshow(img_denoised,[]); title('Denoised (NLM)');

subplot(1,4,3);
imshow(img_dog,[]); title('DoG (should look clean)');

subplot(1,4,4);
imshow(img,[]); hold on;
L_vis = labelmatrix(CC_kept);
overlay = label2rgb(L_vis, 'hsv', 'k');
h = imshow(overlay);
set(h, 'AlphaData', 0.45*(L_vis > 0));
B = bwboundaries(L_vis > 0);
for k = 1:numel(B)
    plot(B{k}(:,2), B{k}(:,1), 'r', 'LineWidth', 1.2);
end
title(sprintf('Detected: %d chromosomes', numel(props_kept)));

%% --- 11. Print stats to help tune filters ---
if numel(props_kept) > 0
    fprintf('\n=== Detection Summary ===\n');
    fprintf('Count : %d\n', numel(props_kept));
    fprintf('Area  : min=%.0f  mean=%.0f  max=%.0f\n', ...
        min(area(keep)), mean(area(keep)), max(area(keep)));
    fprintf('Length: min=%.1f  mean=%.1f  max=%.1f\n', ...
        min(majLen(keep)), mean(majLen(keep)), max(majLen(keep)));
    fprintf('Width : min=%.1f  mean=%.1f  max=%.1f\n', ...
        min(minWid(keep)), mean(minWid(keep)), max(minWid(keep)));

    % Show what was REJECTED to help tune filters
    rejected = props(~keep);
    if numel(rejected) > 0
        r_area = [rejected.Area];
        r_len  = [rejected.MajorAxisLength];
        fprintf('\nRejected %d objects — area range [%.0f–%.0f], length [%.1f–%.1f]\n', ...
            numel(rejected), min(r_area), max(r_area), min(r_len), max(r_len));
    end
end
end
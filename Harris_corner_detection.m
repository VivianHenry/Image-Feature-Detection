clear;
clc;
close all;

%% Reading an image

url1 = "/Users/vivianhenry/Documents/Purdue/Sem 2/Image Processing and Computer Vision/Programming Assignments/3/building.pnm";
url2 = "/Users/vivianhenry/Documents/Purdue/Sem 2/Image Processing and Computer Vision/Programming Assignments/3/hinge.pnm";
url3 = "/Users/vivianhenry/Documents/Purdue/Sem 2/Image Processing and Computer Vision/Programming Assignments/3/hinges.pnm";
url4 = "/Users/vivianhenry/Documents/Purdue/Sem 2/Image Processing and Computer Vision/Programming Assignments/3/keys.pnm";
url5 = "/Users/vivianhenry/Documents/Purdue/Sem 2/Image Processing and Computer Vision/Programming Assignments/3/pillsetc.pnm";
image = imread(url1);
image_original = image(:, :);

%% Attain gradients

%Prewitt mask is used
row_fil = [-1, 0, 1; -1, 0, 1; -1, 0, 1];
colomn_fil = row_fil';

Ix = imfilter(double(image_original), row_fil);
Iy = imfilter(double(image_original), colomn_fil);

%% Gaussian filter to smoothen any noise

filter = fspecial("gaussian", 3, 1);

%% Calculate entries of M matrix

Ix2 = imfilter(double(Ix .^ 2), filter, "same");
Iy2 = imfilter(double(Iy .^ 2), filter, "same");
IxIy = imfilter(double(Ix .* Iy), filter, "same");

%% Harris corner measure

R = ((Ix2 .* Iy2) - (IxIy .^2)) ./ (Ix2 - Iy2);

%% Find local maximum

area_size = 1;
order = (2 * area_size + 1);
local_max = ordfilt2(R, order^2, ones(order));  % Used for 2-D statistical filtering
threshold = 1000000;
harris_points = (R == local_max) & (abs(R) > threshold);  % Returns 1 if satisfied and 0 if not
[rows, cols, val] = find(harris_points);  % Returns index and value of any non-zero entry
count_corners_user = size(val, 1);

%% Plotting

subplot(1, 3, 1);
imshow(image_original); title("Original Image");

subplot(1, 3, 2);
imshow(image_original); hold on;
plot(cols, rows, "ys");
title("User Function");

corners = detectHarrisFeatures(image_original);
count_corners_built_in = size(corners, 1);
subplot(1, 3, 3);
imshow(image_original); hold on;
plot(corners);
title("Built-in Function");

%% Effect of increasing threshold value

Y = [];
index = 0;
for thresh = 2000:1000:300000
    harris_points = (R == local_max) & (abs(R) > thresh);
    [rows, cols, val] = find(harris_points);
    count_corners_user = size(val, 1);
    index = index + 1;
    Y(index) = abs(count_corners_built_in - count_corners_user) / count_corners_built_in;
end

figure;
X = 2000:1000:300000;
plot(X, Y); xlabel("Threshold value"); ylabel("% diff in corner count between built-in and user defined functions")
title("Effect of increasing threshold value");

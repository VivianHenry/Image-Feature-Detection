clear;
clc;
close all;

%% Reading an image

url1 = "/Users/vivianhenry/Documents/Purdue/Sem 2/Image Processing and Computer Vision/Programming Assignments/3/building.pnm";
url2 = "/Users/vivianhenry/Documents/Purdue/Sem 2/Image Processing and Computer Vision/Programming Assignments/3/hinge.pnm";
url3 = "/Users/vivianhenry/Documents/Purdue/Sem 2/Image Processing and Computer Vision/Programming Assignments/3/hinges.pnm";
url4 = "/Users/vivianhenry/Documents/Purdue/Sem 2/Image Processing and Computer Vision/Programming Assignments/3/keys.pnm";
url5 = "/Users/vivianhenry/Documents/Purdue/Sem 2/Image Processing and Computer Vision/Programming Assignments/3/pillsetc.pnm";
image = imread(url4);
image_original = image(:, :);  % Incase the image has 3 dimensions
image_original = medfilt2(image_original, [3, 3]);

%% Canny edge detection

image_original = im2gray(image_original);  % Incase the image is not grayscale
image_original = double (image_original);

% Gaussian smoothening
gauss_filter = fspecial("gaussian", 5, 0.5);
A=conv2(image_original, gauss_filter, 'same');

% Derivative masks for X and Y direction
% Sobel mask is used
dx = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
dy = [-1, -2, -1; 0, 0, 0; 1, 2, 1];

%Convolution of image and derivative masks
Ix = conv2(A, dx, 'same');
Iy = conv2(A, dy, 'same');

% Angle calculation
angle = atan2(Iy, Ix);
angle = angle * (180 / pi);  % Radian to degree conversion

[r, c] = size(A);

%Adjustment for negative directions, making all directions positive
for i = 1:r
    for j = 1:c
        if (angle(i, j) < 0) 
            angle(i,j) = 360 + angle(i, j);
        end
    end
end

modified_angles = zeros(r, c);

%Adjusting directions to nearest 0, 45, 90, or 135 degree
for i = 1:r
    for j = 1:c
        if ((angle(i, j) >= 0 ) && (angle(i, j) < 22.5) || (angle(i, j) >= 157.5) && (angle(i, j) < 202.5) || (angle(i, j) >= 337.5) && (angle(i, j) <= 360))
            modified_angles(i, j) = 0;
        elseif ((angle(i, j) >= 22.5) && (angle(i, j) < 67.5) || (angle(i, j) >= 202.5) && (angle(i, j) < 247.5))
            modified_angles(i, j) = 45;
        elseif ((angle(i, j) >= 67.5 && angle(i, j) < 112.5) || (angle(i, j) >= 247.5 && angle(i, j) < 292.5))
            modified_angles(i, j) = 90;
        elseif ((angle(i, j) >= 112.5 && angle(i, j) < 157.5) || (angle(i, j) >= 292.5 && angle(i, j) < 337.5))
            modified_angles(i, j) = 135;
        end
    end
end

%Calculate magnitude
rms_of_der = (Ix.^2) + (Iy.^2);
rms_of_der = sqrt(rms_of_der);

%Non-Maximum Supression
non_max_supp = zeros (r, c);

% Loop returns a binary matrix with 1 in the position of local maximums
for i = 2:(r - 1)
    for j = 2 :(c - 1)
        if modified_angles(i, j) == 0
            non_max_supp(i, j) = (rms_of_der(i, j) == max([rms_of_der(i, j), rms_of_der(i, j + 1), rms_of_der(i, j - 1)]));
        elseif (modified_angles(i, j) == 45)
            non_max_supp(i, j) = (rms_of_der(i, j) == max([rms_of_der(i, j), rms_of_der(i + 1, j - 1), rms_of_der(i - 1, j + 1)]));
        elseif (modified_angles(i, j) == 90)
            non_max_supp(i, j) = (rms_of_der(i, j) == max([rms_of_der(i, j), rms_of_der(i + 1, j), rms_of_der(i - 1, j)]));
        elseif (modified_angles(i, j) == 135)
            non_max_supp(i, j) = (rms_of_der(i, j) == max([rms_of_der(i, j), rms_of_der(i + 1, j + 1), rms_of_der(i - 1, j - 1)]));
        end
    end
end

non_max_supp = non_max_supp.*rms_of_der;  % Inserts RMS value at local maximums

%Hysteresis Thresholding

threshold_low_frac = 0.075;
threshold_high_frac = 0.175;
threshold_low = threshold_low_frac * max(max(non_max_supp));
threshold_high = threshold_high_frac * max(max(non_max_supp));

image_final = zeros (r, c);

for i = 1  : r
    for j = 1 : c
        if (non_max_supp(i, j) < threshold_low)
            image_final(i, j) = 0;
        elseif (non_max_supp(i, j) > threshold_high)
            image_final(i, j) = 1;
        %Using 8-connected components
        elseif (non_max_supp(i + 1, j) > threshold_high || non_max_supp(i - 1, j) > threshold_high || non_max_supp(i, j + 1) > threshold_high || non_max_supp(i, j - 1) > threshold_high || non_max_supp(i - 1, j - 1) > threshold_high || non_max_supp(i - 1, j + 1) > threshold_high || non_max_supp(i + 1, j + 1) > threshold_high || non_max_supp(i + 1, j - 1) > threshold_high)
            image_final(i,j) = 1;
        end
    end
end

%% Built in function

image_built_in = edge(image_original,'canny');

%% Plotting

figure;
subplot(2, 2, 1);
imshow(uint8(image_original)); title("Original Image");

subplot(2, 2, 2);
imshow(non_max_supp); title("Image after non-max suppression");

subplot(2, 2, 3);
imshow(uint8(image_final.*255)); title("Final result of Canny method (User)");

subplot(2, 2, 4);
imshow(image_built_in); title("Final result of Canny method (Built-in)");

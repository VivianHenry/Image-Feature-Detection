% Input is a binary image, such as that attained after Canny
% Need to run Canny edge detection before running this program
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

%% Hough transform (User; by using the edge image)

image_initial_from_Canny = image_final;

[y,x] = find(image_initial_from_Canny);
[sy,sx] = size(image_initial_from_Canny);
total_pixels = length(x);
maxrho = round(sqrt(sx^2 + sy^2));

hough_matrix_from_edge = zeros(2 * maxrho, 180);

for i = 1:total_pixels
    j = 1;
    for theta = (-pi / 2):(pi / 180):((pi / 2) - (pi / 180))
        rho = round(x(i) .* cos(theta) + y(i) .* sin(theta));
        hough_matrix_from_edge(rho + maxrho, j) = hough_matrix_from_edge(rho + maxrho, j) + 1;
        j = j + 1;
    end
end

%% Hough transform (User; by taking an original image and converting it into binary)

image_initial = imbinarize(image_original);

[y,x]=find(image_initial);
[sy,sx]=size(image_initial);
total_pixels = length(x);
maxrho = round(sqrt(sx^2 + sy^2));

hough_matrix_from_binary = zeros(2*maxrho,180);

for i = 1:total_pixels
    j = 1;
    for theta = (-pi / 2):(pi / 180):((pi / 2) - (pi / 180))
        rho = round(x(i) .* cos(theta) + y(i) .* sin(theta));
        hough_matrix_from_binary(rho + maxrho, j) = hough_matrix_from_binary(rho + maxrho, j) + 1;
        j = j + 1;
    end
end

%% Built in function

image_built_in = edge(image_original,'canny');
[H, T, R] = hough(image_built_in, "RhoResolution", 0.5, "Theta", -90:0.5:89);

%% Plotting

subplot(3, 2, 1);
imshow(image_original); title("Original image");

subplot(3, 2, 2);
imshow(image_initial_from_Canny); title("Edge image from Canny");

subplot(3, 2, 3);
theta = rad2deg(-pi/2:pi/180:pi/2-pi/180);
rho = -maxrho:maxrho-1;
imshow(uint8(hough_matrix_from_edge),[],'xdata',theta,'ydata',rho); title("Hough transform from edge image");
xlabel("Theta"); ylabel("Rho");
axis on, axis normal; hold on; colormap(gca, hot);

subplot(3, 2, 4);
imshow(image_initial); title("Binary image");

subplot(3, 2, 5);
theta = rad2deg(-pi/2:pi/180:pi/2-pi/180);
rho = -maxrho:maxrho-1;
imshow(uint8(hough_matrix_from_binary),[],'xdata',theta,'ydata',rho); title("Hough transfrom from binary image");
xlabel("Theta"); ylabel("Rho");
axis on, axis normal; hold on; colormap(gca, hot);

subplot(3, 2, 6);
imshow(imadjust(rescale(H)), "XData", T, "YData", R, "InitialMagnification", "fit"); title("Hough transform (built-in)");
xlabel("Theta"); ylabel("Rho");
axis on; axis normal; hold on; colormap(gca, hot);

% Need to run Canny edge detection and Hough transform before running this program
clc;
close all;

%% Reading an image

url1 = "/Users/vivianhenry/Documents/Purdue/Sem 2/Image Processing and Computer Vision/Programming Assignments/3/building.pnm";
url2 = "/Users/vivianhenry/Documents/Purdue/Sem 2/Image Processing and Computer Vision/Programming Assignments/3/hinge.pnm";
url3 = "/Users/vivianhenry/Documents/Purdue/Sem 2/Image Processing and Computer Vision/Programming Assignments/3/hinges.pnm";
url4 = "/Users/vivianhenry/Documents/Purdue/Sem 2/Image Processing and Computer Vision/Programming Assignments/3/keys.pnm";
url5 = "/Users/vivianhenry/Documents/Purdue/Sem 2/Image Processing and Computer Vision/Programming Assignments/3/pillsetc.pnm";
image = imread(url3);
image_original = image(:, :);  % Incase the image has 3 dimensions

%% Line detection and fitting (User)

% Hough bin matrix

edge_image = image_final;
rhoStep = 1;  % Quantization of rho
theta = -90:89;  % Quantization of theta

max_length = sqrt((size(edge_image, 1) - 1) ^ 2 + (size(edge_image, 2) - 1) ^ 2);
rho = 2 * (ceil(max_length / rhoStep)) + 1;
diagonal = rhoStep * ceil(max_length / rhoStep); % rho ranges from -diagonal to diagonal
numTheta = length(theta);

Harris_matrix = zeros(rho, numTheta);
rho = -diagonal : diagonal;

for i = 1 : size(edge_image, 1)
    for j = 1 : size(edge_image, 2)
        if edge_image(i, j)
            for k = 1 : numTheta
                temp = j * cos(theta(k) * pi / 180) + i * sin(theta(k) * pi / 180);
                row_index = round((temp + diagonal) / rhoStep) + 1;
                Harris_matrix(row_index, k) = Harris_matrix(row_index, k) + 1;                   
            end
        end            
    end
end

%Hough peaks

numpeaks = 20;  % can be varied at user discretion
threshold = 0.50 * max(Harris_matrix(:));  % 50% of max value
n_hood_size = floor(size(Harris_matrix) / 100.0) * 2 + 1;

peaks = zeros(numpeaks, 2);
num = 0;

while num < numpeaks
    max_value = max(Harris_matrix(:));

    if max_value >= threshold
        num = num + 1;
        [r, c] = find(Harris_matrix == max_value);
        peaks(num, :) = [r(1), c(1)];
        r_start = max(1, r - (n_hood_size(1) - 1) / 2);
        r_end = min(size(Harris_matrix,1), r + (n_hood_size(1) - 1) / 2);
        c_start = max(1, c - (n_hood_size(2) - 1) / 2);
        c_end = min(size(Harris_matrix,2), c + (n_hood_size(2) - 1) / 2);
        for i = r_start : r_end
            for j = c_start : c_end
                    Harris_matrix(i,j) = 0;
            end
        end
    else
        break;          
    end
end

peaks = peaks(1:num, :);

%% Line detection and fitting (Built-in)

BW = edge(image_original,'canny');
[H,T,R] = hough(BW);
P  = houghpeaks(H, numpeaks,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(BW,T,R,P,'FillGap',5,'MinLength',7);

%% Plotting

subplot(2, 2, 1);  % Original image
imshow(image_original); title("Original image");

subplot(2, 2, 2);  % Edge image
imshow(edge_image); title("Edge image");

subplot(2, 2, 3);  % Lines from user code
imshow(image_original); title("Line detection and fitting (User)");
hold on;
for i = 1:size(peaks, 1)
   rho_i = rho(peaks(i, 1));
   theta_i = theta(peaks(i, 2)) * (pi / 180);
   if theta_i == 0
       x1 = rho_i;
       x2 = rho_i;
       if rho_i > 0
           y1 = 1;
           y2 = size(image_original,1);
           plot([x1, x2], [y1, y2], 'r', 'LineWidth', 2); 
       end           
   else
       x1 = 1;
       x2 = size(image_original, 2);
       y1 = (rho_i - x1 * cos(theta_i)) / sin(theta_i);
       y2 = (rho_i - x2 * cos(theta_i)) / sin(theta_i);
       plot([x1, x2], [y1, y2], 'r', 'LineWidth', 2); 
   end         
end

subplot(2, 2, 4);  % Lines from built-in function
imshow(image_original); title("Line detection and fitting (Built-in)");
hold on;
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end



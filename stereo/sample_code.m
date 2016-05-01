%%
clc;
clear;
%% load images and match files for the first example
%%

I1 = imread('house1.jpg');
I2 = imread('house2.jpg');
matches = load('house_matches.txt'); 

I1 = imread('library1.jpg');
I2 = imread('library2.jpg');
matches = load('library_matches.txt'); 

% this is a N x 4 file where the first two numbers of each row
% are coordinates of corners in the first image and the last two
% are coordinates of corresponding corners in the second image: 
% matches(i,1:2) is a point in the first image
% matches(i,3:4) is a corresponding point in the second image

N = size(matches,1);

img_left = im2double(rgb2gray(imread('house1.jpg')));
img_right = im2double(rgb2gray(imread('house2.jpg')));
img_left = im2double(rgb2gray(imread('library1.jpg')));
img_right = im2double(rgb2gray(imread('library2.jpg')));

%% HARRIS DETECTOR
sigma = 3;
thresh = 0.001;
radius = 1;
[cim1, r1, c1] = harris(img_left, sigma, thresh, radius, 0);
p1 = [r1' ; c1'];

[cim2, r2, c2] = harris(img_right, sigma, thresh, radius, 0);
p2 = [r2' ; c2'];


%% NORMALIZED CORRELATION
corr_mat = correlatiomatrix(img_left, p1, img_right, p2, 55);
[numrows,numcols] = size(corr_mat);

[~, col] = max(corr_mat,[],2);
[~, row] = max(corr_mat,[],1);    

% Find consistent matches
n = 1:numrows;
ind = find(row(col(n)) == n);
p1ind = ind;
p2ind = col(ind);

% Extract matched points
m1 = p1(:,p1ind);  
m2 = p2(:,p2ind);    

matches_putative = [m1(2,:)' m1(1,:)' m2(2,:)' m2(1,:)'];

figure;
imshow([I1 I2]); hold on;
plot(matches_putative(:,1), matches_putative(:,2), '+r');
plot(matches_putative(:,3)+size(I1,2), matches_putative(:,4), '+r');
%line([matches_putative(:,1) matches_putative(:,3) + size(I1,2)]', matches_putative(:,[2 4])', 'Color', 'r');

%% RANSAC
x1 = [matches_putative(:,1) matches_putative(:,2)];
x2 = [matches_putative(:,3) matches_putative(:,4)];


threshold = 400;
[avg_residual,inliers,pt1,pt2,closest_pt] = ransac(x1,x2,8,threshold,1000);
disp(avg_residual);
% display points and segments of corresponding epipolar lines
figure;
imshow(I2); hold on;
plot(matches_putative(:,3), matches_putative(:,4), '+r');
line([matches_putative(:,3) closest_pt(:,1)]', [matches_putative(:,4) closest_pt(:,2)]', 'Color', 'r');
line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');


%%
%% display two images side-by-side with matches
%% this code is to help you visualize the matches, you don't need
%% to use it to produce the results for the assignment
%%
figure;
imshow([I1 I2]); hold on;
plot(matches(:,1), matches(:,2), '+r');
plot(matches(:,3)+size(I1,2), matches(:,4), '+r');
line([matches(:,1) matches(:,3) + size(I1,2)]', matches(:,[2 4])', 'Color', 'r');


%%
%% display second image with epipolar lines reprojected 
%% from the first image
%%

% first, fit fundamental matrix to the matches
F_unnorm = fit_fundamental(matches,0); % this is a function that you should write
L = (F_unnorm * [matches(:,1:2) ones(N,1)]')'; % transform points from 
% the first image to get epipolar lines in the second image

% find points on epipolar lines L closest to matches(:,3:4)
L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
pt_line_dist = sum(L .* [matches(:,3:4) ones(N,1)],2);
closest_pt = matches(:,3:4) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

% Find mean squared distance
mean_sd1 = sum((matches(:,3) - closest_pt(:,1)).^2 + ...
    (matches(:,4) - closest_pt(:,2)).^2) / size(matches,1);
disp(mean_sd1);

% find endpoints of segment on epipolar line (for display purposes)
pt1 = closest_pt - [L(:,2) -L(:,1)] * 10; % offset from the closest point is 10 pixels
pt2 = closest_pt + [L(:,2) -L(:,1)] * 10;

% display points and segments of corresponding epipolar lines
figure;
subplot(121);
imshow(I2); hold on;
plot(matches(:,3), matches(:,4), '+r');
line([matches(:,3) closest_pt(:,1)]', [matches(:,4) closest_pt(:,2)]', 'Color', 'r');
line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');


F_norm = fit_fundamental(matches,1); % this is a function that you should write
L = (F_norm * [matches(:,1:2) ones(N,1)]')'; % transform points from 
% the first image to get epipolar lines in the second image

% find points on epipolar lines L closest to matches(:,3:4)
L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
pt_line_dist = sum(L .* [matches(:,3:4) ones(N,1)],2);
closest_pt = matches(:,3:4) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

% Find mean squared distance
mean_sd2 = sum((matches(:,3) - closest_pt(:,1)).^2 + ...
    (matches(:,4) - closest_pt(:,2)).^2) / size(matches,1);
disp(mean_sd2);

% find endpoints of segment on epipolar line (for display purposes)
pt1 = closest_pt - [L(:,2) -L(:,1)] * 10; % offset from the closest point is 10 pixels
pt2 = closest_pt + [L(:,2) -L(:,1)] * 10;

subplot(122);
imshow(I2); hold on;
plot(matches(:,3), matches(:,4), '+r');
line([matches(:,3) closest_pt(:,1)]', [matches(:,4) closest_pt(:,2)]', 'Color', 'r');
line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');

%%

P1 = load('house1_camera.txt');
P2 = load('house2_camera.txt');

I1 = imread('house1.jpg');
I2 = imread('house2.jpg');
matches = load('house_matches.txt'); 
%{
P1 = load('library1_camera.txt');
P2 = load('library2_camera.txt');
 
I1 = imread('library1.jpg');
I2 = imread('library2.jpg');
matches = load('library_matches.txt'); 
%}

[U,S,V] = svd(P1);
center1 = V(:,end);
[U,S,V] = svd(P2);
center2 = V(:,end);
for i = 1:4
	center1(i) = center1(i)/center1(4);
    center2(i) = center2(i)/center2(4);
end

x1 = [matches(:,1) matches(:,2)];
x2 = [matches(:,3) matches(:,4)];

[Xw,residual] = triang(P1, x1, P2, x2);
disp(residual);    

%animateplot(Xw', center1, center2);

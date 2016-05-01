clc;
clear;

%% Read Camera Intrinsic Params
load('cam_param');
%K = cameraParams.IntrinsicMatrix;   % Extract Intrinsic parameters
%cxy = K(3,1:2);
%K(1:2,3) = cxy';
%K(3,1:2) = [0 0];


%% Read Image
img1 = im2double(rgb2gray(imread('image1.jpg')));
img2 = im2double(rgb2gray(imread('image2.jpg')));
img3 = im2double(rgb2gray(imread('3.jpg')));
rows = min([size(img1,1),size(img2,1),size(img3,1)]);
cols = min([size(img1,2),size(img2,2),size(img3,2)]);
img1 = imresize(img1,[rows,cols]);
img2 = imresize(img2,[rows,cols]);
img3 = imresize(img3,[rows,cols]);


%% HARRIS DETECTOR
p1 = detectMinEigenFeatures(img1, 'MinQuality', 0.01);
p2 = detectMinEigenFeatures(img2, 'MinQuality', 0.01);
[features1,valid_points1] = extractFeatures(img1,p1);
[features2,valid_points2] = extractFeatures(img2,p2);
indexPairs = matchFeatures(features1,features2);
matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);
figure; showMatchedFeatures(img1,img2,matchedPoints1,matchedPoints2);


%{
%% NORMALIZED CORRELATION
corr_mat = correlatiomatrix(img1, p1, img2, p2, 3);
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

matches = [m1(2,:)' m1(1,:)' m2(2,:)' m2(1,:)'];
m1 = m1';
m2 = m2';

%{
%% Convert to Homogenous coord and normalize
p1 = [p1 ones(size(p1,1),1)]';
p2 = [p2 ones(size(p2,1),1)]';
p3 = [p3 ones(size(p3,1),1)]';
p1_norm = inv(K) * p1;
p2_norm = inv(K) * p2;
p3_norm = inv(K) * p3;

p1_norm = [p1_norm(2,:) ;p1_norm(1,:) ;p1_norm(3,:)];
p2_norm = [p2_norm(2,:) ;p2_norm(1,:) ;p2_norm(3,:)];
p3_norm = [p3_norm(2,:) ;p3_norm(1,:) ;p3_norm(3,:)];
%}

%% Putative matches
figure;
imshow([img1 img2]); hold on;
plot(matches(:,1), matches(:,2), '+r');
plot(matches(:,3)+size(img1,2), matches(:,4), '+r');
line([matches(:,1) matches(:,3) + size(img1,2)]', matches(:,[2 4])', 'Color', 'r');
%}
%%
% Estimate the fundamental matrix
[fRANSAC, inliers] = estimateFundamentalMatrix(matchedPoints1, ...
    matchedPoints2,'Method','RANSAC','NumTrials',2000,...
    'DistanceThreshold',1e-2);

% Find epipolar inliers
inlierPoints1 = matchedPoints1(inliers, :);
inlierPoints2 = matchedPoints2(inliers, :);

% Display inlier matches
figure
showMatchedFeatures(img1, img2, inlierPoints1, inlierPoints2);
title('Epipolar Inliers');

%%
[R, t] = cameraPose(fRANSAC, cameraParams, inlierPoints1, inlierPoints2);

camMatrix1 = cameraMatrix(cameraParams, eye(3), [0 0 0]);
camMatrix2 = cameraMatrix(cameraParams, R', -t*R');

p1 = detectMinEigenFeatures(img1, 'MinQuality', 0.008);
p2 = detectMinEigenFeatures(img2, 'MinQuality', 0.008);
[features1,valid_points1] = extractFeatures(img1,p1);
[features2,valid_points2] = extractFeatures(img2,p2);
indexPairs = matchFeatures(features1,features2);
matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);

% Visualize the camera locations and orientations
cameraSize = 0.3;
figure
plotCamera('Size', cameraSize, 'Color', 'r', 'Label', '1', 'Opacity', 0);
hold on
grid on
plotCamera('Location', t, 'Orientation', R, 'Size', cameraSize, ...
    'Color', 'b', 'Label', '2', 'Opacity', 0);


% Compute the 3-D points
points3D = triangulate(matchedPoints1, matchedPoints2, camMatrix1, camMatrix2);
%figure;
%plot3(points3D(:,1),points3D(:,2),points3D(:,3),'+');
%axis equal;

ptCloud = pointCloud(points3D);
% Visualize the point cloud
pcshow(ptCloud, 'VerticalAxis', 'y', 'VerticalAxisDir', 'down', ...
    'MarkerSize', 200);
camorbit(0, -30);
camzoom(1.5);

%% left: left image matching points
%% right: right image matching points
%% num_samples: min number of samples to choose from
%% threshold: max distance for transformation
%% num_cycles: Loop N times

function [avg_residual,inliers,pt1,pt2,closest_pt] = ransac(left,right,num_samples,threshold,num_cycles)

n = 1:size(left,1);
transformed_pts = zeros(size(left,1),2);
num_inliers = 0;

for i = 1:num_cycles
    ind = datasample(n,num_samples,'Replace',false)';
    left_sample = left(ind,:);
    right_sample = right(ind,:);
    matches = [left_sample right_sample];
    F = fit_fundamental(matches,1);
    N = size(right,1);
    L = (F * [left ones(N,1)]')'; % transform points from 
    % the first image to get epipolar lines in the second image
    % find points on epipolar lines L closest to matches(:,3:4)
    L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
    pt_line_dist = sum(L .* [right ones(N,1)],2);
    closest_pt = right - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

    dist = sum((right - closest_pt).^2,2);
    inliers_ind = find(dist < threshold);
    if numel(inliers_ind) > num_inliers
        num_inliers = numel(inliers_ind);
        inliers = inliers_ind;
        avg_residual = sum(dist) / size(left,1);
    end
end

left_sample = left(inliers,:);
right_sample = right(inliers,:);
matches = [left_sample right_sample];
F = fit_fundamental(matches,1);
N = size(right,1);
L = (F * [left ones(N,1)]')'; % transform points from 
% the first image to get epipolar lines in the second image
% find points on epipolar lines L closest to matches(:,3:4)
L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
pt_line_dist = sum(L .* [right ones(N,1)],2);
closest_pt = right - L(:,1:2) .* repmat(pt_line_dist, 1, 2);
dist = sum((right - closest_pt).^2,2);
inliers_ind = find(dist < threshold);
num_inliers = numel(inliers_ind);
inliers = inliers_ind;
avg_residual = sum((right(:,1) - closest_pt(:,1)).^2 + ...
    (right(:,2) - closest_pt(:,2)).^2) / size(right,1);
% find endpoints of segment on epipolar line (for display purposes)
pt1 = closest_pt - [L(:,2) -L(:,1)] * 10; % offset from the closest point is 10 pixels
pt2 = closest_pt + [L(:,2) -L(:,1)] * 10;


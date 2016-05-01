% Copyright (c) 2004-2009 Peter Kovesi
% Adapted and Modified copy
% Modified different implementation of finding normalized correlation

function corr_mat = correlatiomatrix(im1, p1, im2, p2, w)

% Initialize correlation matrix values to -infinty
corr_mat = -ones(size(p1,2),size(p2,2))*Inf;

[im1rows, im1cols] = size(im1);
[im2rows, im2cols] = size(im2);    

r = (w-1)/2;   

% Remove boundary indices
n1ind = find(p1(1,:)>r & p1(1,:)<im1rows+1-r & ...
    p1(2,:)>r & p1(2,:)<im1cols+1-r);

n2ind = find(p2(1,:)>r & p2(1,:)<im2rows+1-r & ...
    p2(2,:)>r & p2(2,:)<im2cols+1-r);    

for n1 = n1ind            
    % Generate window in 1st image   	
    w1 = im1(p1(1,n1)-r:p1(1,n1)+r, p1(2,n1)-r:p1(2,n1)+r);
    w1 = reshape(w1,1,size(w1,1)*size(w1,2));

    % Calculate noralised correlation measure
    for n2 = n2ind 
        % Generate window in 2nd image
        w2 = im2(p2(1,n2)-r:p2(1,n2)+r, p2(2,n2)-r:p2(2,n2)+r);
        w2 = reshape(w2,1,size(w2,1)*size(w2,2));
        corr_mat(n1,n2) = sum(w1.*w2)./sqrt(sum(w1.^2)*sum(w2.^2));
    end
end
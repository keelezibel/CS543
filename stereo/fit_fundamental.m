%% index = 0, Un-normalized
%% index = 1, Normalized

function F = fit_fundamental(matches,index)

if index == 0
    x = matches(:,1);
    x_prime = matches(:,3);
    y = matches(:,2);
    y_prime = matches(:,4);

    F = [x_prime.*x x_prime.*y x_prime y_prime.*x y_prime.*y ...
        y_prime x y ones(numel(x),1)];

    [U,S,V]=svd(F); 
    F = V(:,end);               
    F = reshape(F,3,3)';        % Reshape to 3x3
    [U,S,V] = svd(F);
    S(end) = 0;                 % Setting smallest to 0
    F = U * S * V';
else
    %% Normalize
    num_matches = size(matches,1);
    average_pixel = sum(matches(:,1:2),1)./num_matches;
    match_norm = matches(:,1:2) - repmat(average_pixel,num_matches,1);    % Zero mean
    mean_sd = sum(sum(match_norm.^2,2))/num_matches;   % mean squared distance
    scaling1 = (2./mean_sd);      % Scale to 2 pixels
    % Line 25: https://www.youtube.com/watch?v=NKxXGsZdDp8
    T1 = [(scaling1)*eye(2), (-(scaling1)*average_pixel)'; 0 0 1];
    match_norm = [T1 * [matches(:,1:2) ones(num_matches,1)]']';

    average_pixel = sum(matches(:,3:4),1)./num_matches;
    match_norm2 = matches(:,3:4) - repmat(average_pixel,num_matches,1);    % Zero mean
    mean_sd = sum(sum(match_norm2.^2,2))/num_matches;   % mean squared distance
    scaling2 = (2./mean_sd);      % Scale to 2 pixels
    % Line 33: https://www.youtube.com/watch?v=NKxXGsZdDp8
    T2 = [(scaling2)*eye(2), (-(scaling2)*average_pixel)'; 0 0 1];
    match_norm2 = [T2 * [matches(:,3:4) ones(num_matches,1)]']';

    x = match_norm(:,1);
    x_prime = match_norm2(:,1);
    y = match_norm(:,2);
    y_prime = match_norm2(:,2);

    F = [x_prime.*x x_prime.*y x_prime y_prime.*x y_prime.*y ...
        y_prime x y ones(numel(x),1)];

    [U,S,V]=svd(F); 
    F = V(:,end);               
    F = reshape(F,3,3)';        % Reshape to 3x3
    [U,S,V] = svd(F);
    S(end) = 0;                 % Setting smallest to 0
    F_norm = U * S * V';
    F=T2'*F_norm*T1;
end

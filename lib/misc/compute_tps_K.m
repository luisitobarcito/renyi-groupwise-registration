function [K] = compute_tps_K(x, y)
% [K] = compute_tps_K(x, y)
%%
% Purpose: Generate K (TPS kernel) matrix   || x - y ||^2 * log || x - y ||
% x = n x d
% y = m x d

if(nargin ~= 2)
    disp('# ERROR #: compute_dist requires two inputs !');
    return;
end

[n, d] = size(x);
[m, d2] = size(y);
if(d ~= d2)
    disp('# ERROR #: dimensions of two variables dont match !');
    return;
end

% compute || x - y ||^2
K = compute_dist(x, y);

% 2D: K = || x - y ||^2 * log || x - y ||
% 3D: K = || x - y ||
if(d == 2)
    % replace points close to 0 with 1 to avoid singularity (log(0) is undefined)
    mask = K < 1e-10;
    % compute K * log( sqrt(K) ) and replace mask samples again with 0
    K = K .* ( (1/2 * log(K + mask)) .* ~mask );
else
    K = sqrt(K);
end

end

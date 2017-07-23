function [K] = compute_dist(x, y)
% [K] = compute_dist(x, y)
%%
% Purpose: compute   || x - y ||^2
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
K = bsxfun(@minus, sum(x.*x, 2), 2*x*y');
K = bsxfun(@plus, K, sum(y'.*y', 1)); 
%diag(x*x')*ones(1,m) - 2*x*y' + ones(n,1)*diag(y*y')';

end

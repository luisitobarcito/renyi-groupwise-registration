function [K] = compute_gauss_K(x, y, sig, normalized)
% [K] = compute_gauss_K(x, y, sig)
%%
% Purpose: Generate K (Gaussian kernel) matrix   exp( -|| x - y ||^2 / 2*sig^2 )
% x = n x d
% y = m x d

if(nargin > 4 || nargin < 3)
    disp('# ERROR #: compute_dist requires three or four inputs !');
    return;
end

if nargin == 3
   normalized = true;
end

[n, d] = size(x);
[m, d2] = size(y);
if(d ~= d2)
    disp('# ERROR #: dimensions of two variables dont match !');
    return;
end

% compute || x - y ||^2
K = compute_dist(x, y);
% compute Gaussian kernel matrix
if normalized
   K = (1/sqrt(2*pi*sig^2)) * exp( -K / (2*sig^2) );
else 
   K = exp( -K / (2*sig^2) );
end
end

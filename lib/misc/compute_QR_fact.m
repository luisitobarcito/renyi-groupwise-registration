function [Q1, Q2, R] = compute_QR_fact(x)
% [Q1, Q2, R] = compute_QR_fact(x)
%%
% Purpose: Generete QR decomposition
%input
% x = n x d
%output
% Q1 = n x d
% Q2 = n x (n-d)
% R  = d x d

if(nargin ~= 1)
    disp('# ERROR #: compute_dist requires one input !');
    return;
end

[n, d] = size(x);

[Q, R] = qr(x);
Q1 = Q(:, 1:d);
Q2 = Q(:, d+1:n);
R  = R(1:d, 1:d);

end

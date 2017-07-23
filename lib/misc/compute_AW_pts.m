function [A, W] = compute_AW_pts(Q1, Q2, R, K, y, lmbd, expVal, sig)
% [A, W] = compute_AW_pts(Q1, Q2, R, K, y, lmbd, expVal, sig)
%%
% Purpose: Calculate A and W transformation matrices
%input
% Q1 = n x d
% Q2 = n x (n-d)
% R  = d x d
% K  = n x n
% y  = n x d
% lmbd = scalar
%output
% A = d x d
% W = n x d

if(nargin ~= 8)
    disp('# ERROR #: compute_dist requires eight inputs !');
    return;
end

[n, d] = size(Q1);
[m, d2] = size(y);
if( (d ~= d2) || (n ~= m) )
    disp('# ERROR #: dimensions of the variables dont match !');
    return;
end

% %     W = Q2 / (Q2'*(K + lmbd*eye(n))*Q2) * Q2'*y;
% %     A = R \ Q1' * (y - (K + lmbd*eye(n))*W);

E = expVal;
% %     De = diag(ones(1,m)*E);
De = diag(E*ones(m,1));
M = K*De + lmbd*eye(n);

W = Q2 / (Q2'*M*Q2) * Q2'*E*y;
A = (De(1:d,1:d)*R) \ Q1' * (E*y - M*W);

end

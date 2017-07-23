function [tS] = sum_over_dims(G, dims)
% [tS] = sum_over_dims(G, dims)
%%
% Purpose: compute   sum(G, [dimensions]) over several dimensions
% G = Nx x Ny x Nz x ....
% dims = [vector of dimensions]
if(nargin ~= 2)
    disp('# ERROR #: compute_dist requires two inputs !');
    return;
end

tS = G;
for d = 1:length(dims)
    tS = sum(tS, dims(d));
end
tS = squeeze(tS);

end

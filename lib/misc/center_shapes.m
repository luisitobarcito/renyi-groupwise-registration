function Xcent = center_shapes(X)
% Xcent = center_shapes(X)
% Center all shapes around the mean of the union of all the points in
% the set of shapes.
tXmean = zeros(1, size(X{1}, 2));
totalX = 0;
for iSh = 1:length(X)
    tXmean = tXmean + size(X{iSh},1)*mean(X{iSh});
    totalX = totalX + size(X{iSh},1);
end
tXmean = tXmean/totalX;
for iSh = 1:length(X)
    Xcent{iSh} = bsxfun(@minus, X{iSh}, tXmean);
end

end 

function overlay_shapes(X, varargin)
% overlay_shapes(X, markers)

if (nargin > 3)
  error('Too many input parameters');
end
if isempty(varargin)
   markers = {'ro', 'bo', 'go', 'mo', 'co', 'yo', 'ko' };
else
   markers = varargin{1};
end
colors = 'rbgmcyk';
hold on
for iSh = 1: length(X)
  marker_type = markers{iSh};
  marker_facecolor = intersect(marker_type, colors);
  plot(X{iSh}(:, 1), X{iSh}(:, 2), marker_type, 'MarkerSize', 5, 'LineWidth', 1, 'MarkerFaceColor', marker_facecolor);
  %plot(X{iSh}(:, 1), X{iSh}(:, 2), marker_type, 'MarkerSize', 10, 'MarkerFaceColor', marker_facecolor);
end
end

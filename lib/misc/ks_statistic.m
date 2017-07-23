function [ KSstatistic] = ks_statistic(x1, x2)

% ks_statistic - FUNCTION Two-sample Two-diensional Kolmogorov-Smirnov statistic
%
% Usage:[KSstatistic] = ks_statistic(x1, x2)
%
% Kolmogorov-Smirnov statistic is used to
% determine whether two sets of data arise from the same or different
% distributions.  The null hypothesis is that both data sets were drawn
% from the same continuous distribution.
% 
% The algorithm in this function is taken from Peacock [1].
% 
% 'x1' is an [Nx2] matrix, each row containing a two-dimensional sample.
% 'x2' is an [Mx2] matrix, each row likewise containing a two-dimensional
% sample. 
% 'KSstatistic' is the raw value for the test statistic ('D' in [1]).
% 
% References: [1] J. A. Peacock, "Two-dimensional goodness-of-fit testing
%  in astronomy", Monthly Notices Royal Astronomy Society 202 (1983)
%  615-627.
%    Available from: http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1983MNRAS.202..615P&defaultprint=YES&filetype=.pdf
%

% Adapted from Dylan Muir  kstest_2s_2d
%%


if ((size(x1,2)~=2)||(size(x2,2)~=2))
   error('stats:kstest2:TwoColumnMatrixRequired','The samples X1 and X2 must be two-column matrices.');
end
n1 = size(x1,1);
n2 = size(x2,1);




%%
%
% Calculate F1(x) and F2(x), the empirical (i.e., sample) CDFs.
%

% - A function handle to perform comparisons in all possible directions
fhCounts = @(x, edge)([(x(:, 1) >= edge(1)) & (x(:, 2) >= edge(2))...
   (x(:, 1) <= edge(1)) & (x(:, 2) >= edge(2))...
   (x(:, 1) <= edge(1)) & (x(:, 2) <= edge(2))...
   (x(:, 1) >= edge(1)) & (x(:, 2) <= edge(2))]);

KSstatistic = -inf;

for iX = 1:(n1+n2)
   % - Choose a starting point
   if (iX<=n1)
      edge = x1(iX,:);
   else
      edge = x2(iX-n1,:);
   end
   
   % - Estimate the CDFs for both distributions around this point
   vfCDF1 = sum(fhCounts(x1, edge)) ./ n1;
   vfCDF2 = sum(fhCounts(x2, edge)) ./ n2;
   
   % - Two-tailed test statistic
   vfThisKSTS = abs(vfCDF1 - vfCDF2);
   fKSTS = max(vfThisKSTS);
   
   % - Final test statistic is the maximum absolute difference in CDFs
   if (fKSTS > KSstatistic)
      KSstatistic = fKSTS;
   end
end

end

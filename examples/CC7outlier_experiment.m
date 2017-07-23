% This script compares the following methods for group-wise
% registration on the Corpus Callosum data set where one of the shapes have
% been contaminated with outliers:
%  1- CDF-HC method described in,
%     "Group-wise point-set registration using a novel CDF-based
%     Havrda-Charvat divergence" by T. Chen, B. C. Vemuri, A.
%     Rangarajan, and S. J. Eisenschenk.
%  2- Holder's inequality, 
%  3- Information potential difference, and
%  4- Scaled information potential differnece.
% 
% The last threr methods are described in,
%     "Group-wise point-set regisytration based on Renyi's second
%     order entropy" by L. Sanchez Giraldo, E. Hasanbelliu, M. Rao,
%     and J. C. Principe.    

clear all
close all
clc

%% Add multishape aligment functions to the search path 
disp('Load shape alignment library files.');
run ../lib/init_lib.m;

%% Load the Corpus Callosum data 
load ../data/cdfhc_data2D_CC.mat;
disp('Loaded Corpus Callosum data')
X = CC7;
%% Load the Corpus Callosum with outliers data 
load ../data/cdfhc_data2D_CC_outliers.mat;
disp('Loaded Corpus Callosum data with outliers')
X_outliers = CC7_outliers;



%% Center all shapes with respect to the origin
X = center_shapes(X);
X_outliers = center_shapes(X_outliers);

%% Display set of centered shapes before registration
hfig = figure();
overlay_shapes(X);
title('Shapes before registration')
xlim([-0.4 0.4])
ylim([-0.2 0.15])
set(gca, 'DataAspectRatio', [2 1 1]);
drawnow, pause(0.1)
close(hfig)

%% Run cdfHC
tXcdfhc = HC2Reg_TPS(X, 0, 300);
tXcdfhc = center_shapes(tXcdfhc);
tXcdfhc_outliers = HC2Reg_TPS(X_outliers, 0, 300);
tXcdfhc_outliers = center_shapes(tXcdfhc_outliers);


%% Run Holder's inequality 
options.disp_progress = 0;
options.regularization = 0.01;
options.anealing = 0.98;
options.max_iter = 300;
tXholder = Holder_FixedPoint(X, options);
tXholder_outliers = Holder_FixedPoint(X_outliers, options);


%% Run Information Potential (IP) differnece and normalized version
% Parameter set same for both algorithms
options.disp_progress = 0;
options.regularization = 0.1;
options.init_sigma = 0.1;
options.final_sigma = 0.01;
options.anealing = 0.98;
options.max_iter = 300;

% unnormalized version 
tXunip = IP_Diff_FixedPoint(X, options);
tXunip_outliers = IP_Diff_FixedPoint(X_outliers, options);

% normalized version
tXnmip = IP_Norm_Diff_FixedPoint(X, options);
tXnmip_outliers = IP_Norm_Diff_FixedPoint(X_outliers, options);


% Plot all original shapes and the outliers  

hfig = figure();
markers = {'ro', 'bo', 'go', 'mo', 'co', 'yo', 'ko', 'k*'};
  shapes = {X_outliers{1:6}, X_outliers{end}(1:size(X_outliers{1}, 1),:),X_outliers{end}(size(X_outliers{1}, 1)+1:end,:) }; 
  overlay_shapes(shapes, markers);
  title('Original shapes and outlier points')
  xlim([-0.4 0.4])
  ylim([-0.2 0.15])
  set(gca, 'DataAspectRatio', [2 1 1]);
drawnow, pause(0.1)
close(hfig);


%% Plot all shapes together for each method on the outlier set and
  %% compare with the original solution witohout outliers 

hfig = figure();
markers = {'ro', 'bo', 'go', 'mo', 'co', 'yo', 'ko', 'b+'};
shapes = {tXcdfhc_outliers{:}, tXcdfhc{end}};
overlay_shapes(shapes, markers);
title('cdf-HC')
xlim([-0.4 0.4])
ylim([-0.2 0.15])
set(gca, 'DataAspectRatio', [2 1 1]);
drawnow, pause(0.1)
close(hfig);


hfig = figure();
shapes = {tXholder_outliers{:}, tXholder{end}};
overlay_shapes(shapes, markers);
title('Holder')
xlim([-0.4 0.4])
ylim([-0.2 0.15])
set(gca, 'DataAspectRatio', [2 1 1]);
drawnow, pause(0.1)
close(hfig);

hfig = figure();
shapes = {tXunip_outliers{:}, tXunip{end}};
overlay_shapes(shapes, markers);
title('IP difference')
xlim([-0.4 0.4])
ylim([-0.2 0.15])
set(gca, 'DataAspectRatio', [2 1 1]);
drawnow, pause(0.1)
close(hfig);

hfig = figure();
shapes = {tXnmip_outliers{:}, tXnmip{end}};
overlay_shapes(shapes, markers);
title('Normalized IP difference')
xlim([-0.4 0.4])
ylim([-0.2 0.15])
set(gca, 'DataAspectRatio', [2 1 1]);
drawnow, pause(0.1)
close(hfig);

%% Compute Kolomogorv-Smirnov  statistic
for i = 2:length(tXunip)
  KSstat_cdfhc(i-1) = ks_statistic(tXcdfhc_outliers{1}(:,1:2), tXcdfhc_outliers{i}(:,1:2));
  KSstat_holder(i-1) = ks_statistic(tXholder_outliers{1}(:,1:2), tXholder_outliers{i}(:,1:2));
  KSstat_unip(i-1) = ks_statistic(tXunip_outliers{1}(:,1:2), tXunip_outliers{i}(:,1:2));
  KSstat_nmip(i-1) = ks_statistic(tXnmip_outliers{1}(:,1:2), tXnmip_outliers{i}(:,1:2));
end
fprintf('(+ outlier) KS_statistic cdf-HC %f\n', mean(KSstat_cdfhc));
fprintf('(+ outlier) KS_statistic Holder %f\n', mean(KSstat_holder));
fprintf('(+ outlier) KS_statistic unnormalized IP %f\n', mean(KSstat_unip));
fprintf('(+ outlier) KS_statistic normalized IP %f\n', mean(KSstat_nmip));

%% Compute Kolomogorv-Smirnov  statistic
for i = 2:length(tXunip)
  KSstat_cdfhc(i-1) = ks_statistic(tXcdfhc_outliers{1}(:,1:2), tXcdfhc_outliers{i}(1:63,1:2));
  KSstat_holder(i-1) = ks_statistic(tXholder_outliers{1}(:,1:2), tXholder_outliers{i}(1:63,1:2));
  KSstat_unip(i-1) = ks_statistic(tXunip_outliers{1}(:,1:2), tXunip_outliers{i}(1:63,1:2));
  KSstat_nmip(i-1) = ks_statistic(tXnmip_outliers{1}(:,1:2), tXnmip_outliers{i}(1:63,1:2));
end
fprintf('(- outlier) KS_statistic cdf-HC %f\n', mean(KSstat_cdfhc));
fprintf('(- outlier) KS_statistic Holder %f\n', mean(KSstat_holder));
fprintf('(- outlier) KS_statistic unnormalized IP %f\n', mean(KSstat_unip));
fprintf('(- outlier) KS_statistic normalized IP %f\n', mean(KSstat_nmip));
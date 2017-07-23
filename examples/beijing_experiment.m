% This script compares the following methods for group-wise
% registration on the Corpus Callosum data set:
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
load ../data/cdfhc_data2D_beijing.mat;
disp('Loaded beijing data')
X = beijing;

%% Center all shapes with respect to the origin
X = center_shapes(X);

%% Display set of centered shapes before registration
hfig = figure();
overlay_shapes(X);
title('Shapes before registration')
xlim([-1.2 1.2])
ylim([-2 1.5])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow, pause(0.1)
close(hfig);

%% Run cdfHC
tXcdfhc = HC2Reg_TPS(X, 1, 300);
tXcdfhc = center_shapes(tXcdfhc);

%% Run Holder's inequality 
options.disp_progress = 0;
options.init_sigma = 0.5;
options.final_sigma = 0.05;
options.regularization = 0.01;
options.anealing = 0.99;
options.max_iter = 300;
tic
tXholder = Holder_FixedPoint(X, options);
toc
%% Run Information Potential (IP) differnece and normalized version
% Parameter set same for both algorithms
options.disp_progress = 0;
options.regularization = 0.01;
options.init_sigma = 0.5;
options.final_sigma = 0.01;
options.anealing = 0.99;
options.max_iter = 300;

% unnormalized version 
tic
tXunip = IP_Diff_FixedPoint(X, options);
toc
% normalized version
tic
tXnmip = IP_Norm_Diff_FixedPoint(X, options);
toc


%% Plot all shapes together for each method  

hfig = figure();
overlay_shapes(tXcdfhc);
title('cdf-HC')
xlim([-1.2 1.2])
ylim([-2 1.5])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow;
close(hfig);

hfig = figure();
overlay_shapes(tXholder);
title('Holder')
xlim([-1.2 1.2])
ylim([-2 1.5])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow;
close(hfig);


hfig = figure();
overlay_shapes(tXunip);
title('IP difference')
xlim([-1.2 1.2])
ylim([-2 1.5])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow;
close(hfig);

hfig= figure();
overlay_shapes(tXnmip);
title('Normalized IP difference')
xlim([-1.2 1.2])
ylim([-2 1.5])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow;
close(hfig);

%% Compute Kolomogorv-Smirnov  statistic
for i = 2:length(tXunip)
  KSstat_cdfhc(i-1) = ks_statistic(tXcdfhc{1}(:,1:2), tXcdfhc{i}(:,1:2));
  KSstat_holder(i-1) = ks_statistic(tXholder{1}(:,1:2), tXholder{i}(:,1:2));
  KSstat_unip(i-1) = ks_statistic(tXunip{1}(:,1:2), tXunip{i}(:,1:2));
  KSstat_nmip(i-1) = ks_statistic(tXnmip{1}(:,1:2), tXnmip{i}(:,1:2));
end
fprintf('KS_statistic cdf-HC %f\n', mean(KSstat_cdfhc));
fprintf('KS_statistic Holder %f\n', mean(KSstat_holder));
fprintf('KS_statistic unnormalized IP %f\n', mean(KSstat_unip));
fprintf('KS_statistic normalized IP %f\n', mean(KSstat_nmip));



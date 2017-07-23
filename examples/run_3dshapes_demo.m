clear all
close all
clc

%% Add multishape aligment functions to the search path 
disp('Load shape alignment library files.');
run ../lib/init_lib.m;

%% Load McGill3D data 
load ../data/mcgill3d_humans.mat;
disp('Loaded 3d data')


% Center all shapes with respect to the origin
tXmean = zeros(1, size(X{1}, 2));
totalX = 0;
for iSh = 1:length(X)
    tXmean = tXmean + size(X{iSh},1)*mean(X{iSh});
    totalX = totalX + size(X{iSh},1);
end

scaleX = 1/128;
tXmean = tXmean/totalX;

% Subsample shapes for computation
n_sub = 2000;
n_shapes = 3;
s_idx = [29 25 14];
% s_idx = randperm(length(X));
for iSh = 1:n_shapes
    p_idx = randperm(size(X{s_idx(iSh)},1));
    Xsub{iSh} = bsxfun(@minus, X{s_idx(iSh)}(p_idx(1:n_sub), :), tXmean)*scaleX;
end
X = Xsub;
%% Display set of shapes before registration
colmark = ['r'; 'b'; 'g'; 'm'; 'c'; 'y'; 'k'];
figure
hold on
for iSh = 1: length(X)
    colidx = mod(iSh-1, length(colmark)) + 1;
    plot3(X{iSh}(:, 1), X{iSh}(:, 2), X{iSh}(:, 3), ['o', colmark(colidx)], 'MarkerSize', 3, 'MarkerFaceColor', colmark(colidx));
end
title('Shapes before registration')
% xlim([-0.4 0.4])
% ylim([-0.4 0.4])
set(gca, 'DataAspectRatio', [1 1 1]);

%% Run unnormalized IP potential difference algorithm
options.disp_progress = 1;
options.init_sigma = 0.1;
options.anealing = 0.98;
options.max_iter = 100;
options.nonrigid =  'Gaussian_RBF';
options.regularization = 0.01;

disp('Running Unnormalized IP algorithm');
[tXun, Aun, Wun] = IP_Diff_FixedPoint(X, options);

%%% Display set of shapes after registration
% Center all shapes with respect to the origin
tXmeanun = zeros(1, size(tXun{1}, 2));
totalXun = 0;
for iSh = 1:length(tXun)
    tXmeanun = tXmeanun + size(tXun{iSh},1)*mean(tXun{iSh});
    totalXun = totalXun + size(tXun{iSh},1);
end
tXmeanun = tXmeanun/totalXun;

% Plot all shpes together 

figure, hold on;
for iSh = 1:length(tXun)
    tXun{iSh} = bsxfun(@minus, tXun{iSh}, tXmeanun);
    colidx = mod(iSh-1, length(colmark)) + 1;
    plot3(tXun{iSh}(:,1), tXun{iSh}(:,2), tXun{iSh}(:,3), ['o',colmark(colidx)], 'MarkerSize',3, 'MarkerFaceColor', colmark(colidx));
end
set(gca, 'DataAspectRatio', [1 1 1]); 
% xlim([-0.4 0.4])
% ylim([-0.4 0.4])
title('Shapes after Registration (Unnormalized IP)');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Run normalized IP potential difference algorithm

disp('Running Normalized IP algorithm');
[tXnm, Anm, Wnm] = IP_Norm_Diff_FixedPoint(X, options);

%%% Display set of shapes after registration
% Center all shapes with respect to the origin
tXmeannm = zeros(1, size(tXnm{1}, 2));
totalXnm = 0;
for iSh = 1:length(tXnm)
    tXmeannm = tXmeannm + size(tXnm{iSh},1)*mean(tXnm{iSh});
    totalXnm = totalXnm + size(tXnm{iSh},1);
end
tXmeannm = tXmeannm/totalXnm;

% Plot all shpes together 
figure, hold on;
for iSh = 1:length(tXnm)
    tXnm{iSh} = bsxfun(@minus, tXnm{iSh}, tXmeannm);
    colidx = mod(iSh-1, length(colmark)) + 1;
    plot3(tXnm{iSh}(:,1), tXnm{iSh}(:,2), tXnm{iSh}(:,3), ['o',colmark(colidx)], 'MarkerSize',3, 'MarkerFaceColor', colmark(colidx));
end
set(gca, 'DataAspectRatio', [1 1 1]); 
xlim([-0.4 0.4])
ylim([-0.4 0.4])
title('Shapes after Registration (Normalized IP)');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


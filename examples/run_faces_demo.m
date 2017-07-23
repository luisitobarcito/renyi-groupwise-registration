clear all
close all
clc

%% Add multishape aligment functions to the search path 
disp('Load shape alignment library files.');
run ../lib/init_lib.m;

%% Load the Corous Callosum data 
load ../data/imm_data2D_faces.mat;
disp('Loaded faces data')
X = faces_data(:);
% Center all shapes with respect to the origin
tXmean = zeros(1, size(X{1}, 2));
totalX = 0;
for iSh = 1:length(X)
    tXmean = tXmean + size(X{iSh},1)*mean(X{iSh});
    totalX = totalX + size(X{iSh},1);
end

scaleX = 1/400;
tXmean = tXmean/totalX;

for iSh = 1:length(X)
    X{iSh} = bsxfun(@minus, X{iSh}, tXmean)*scaleX;
end

%% Display set of shapes before registration
colmark = ['r'; 'b'; 'g'; 'm'; 'c'; 'y'; 'k'];
figure
hold on
for iSh = 1: length(X)
    colidx = mod(iSh-1, length(colmark)) + 1;
    plot(X{iSh}(:, 1), X{iSh}(:, 2), ['o', colmark(colidx)], 'MarkerSize', 10, 'MarkerFaceColor', colmark(colidx));
end
title('Shapes before registration')
xlim([-0.4 0.4])
ylim([-0.4 0.4])
set(gca, 'DataAspectRatio', [1 1 1]);

%% Run unnormalized IP potential difference algorithm
options.disp_progress = 10;
options.init_sigma = 0.1;
options.anealing = 0.98;
options.max_iter = 100;

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
    plot(tXun{iSh}(:,1), tXun{iSh}(:,2), ['o',colmark(colidx)], 'MarkerSize',10, 'MarkerFaceColor', colmark(colidx));
end
set(gca, 'DataAspectRatio', [1 1 1]); 
xlim([-0.4 0.4])
ylim([-0.4 0.4])
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
    plot(tXnm{iSh}(:,1), tXnm{iSh}(:,2), ['*',colmark(colidx)], 'MarkerSize',10, 'MarkerFaceColor', colmark(colidx));
end
set(gca, 'DataAspectRatio', [1 1 1]); 
xlim([-0.4 0.4])
ylim([-0.4 0.4])
title('Shapes after Registration (Normalized IP)');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute Kolomogorv-Smirnov  statistic
%% Compute Kolomogorv-Smirnov  statistic
for i = 2:length(tXun)
  KSstat_un(i-1) = ks_statistic(tXun{1}, tXun{i});
  KSstat_nm(i-1) = ks_statistic(tXnm{1}, tXnm{i});
end

fprintf('KS_statistic unnormalized IP %f\n', mean(KSstat_un));
fprintf('KS_statistic normalized IP %f\n', mean(KSstat_nm));

%% Compute warped images
params.origin = tXmean;
params.scale = scaleX;
params.x_min = -0.6;
params.x_max = 0.6;
params.y_min = -0.6;
params.y_max = 0.6;

im_tot = zeros(480,480,3);
params.size = [480,480,3];
n_im = length(X);    
for iImg = 1:n_im
    fprintf('warping %s\n',annotation_files{iImg}(1:end-4));
    im = imread(strcat('../data/face_data/', annotation_files{iImg}(1:end-4)));
    params.centers = X{iImg};
    params.A = Anm{iImg};
    params.W = Wnm{iImg};
    im_w = warp_image(im, params);
    im_tot = im_tot + im_w/n_im; 
end
figure
imshow(im_tot/255)
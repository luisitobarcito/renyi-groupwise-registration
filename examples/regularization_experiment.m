% This script illusterates the performacce of the group-wise
% shape registration based on the normalized information potential
% difference  described in,
%     "Group-wise point-set registration based on Renyi's second
%     order entropy" by L. Sanchez Giraldo, E. Hasanbelliu, M. Rao,
%     and J. C. Principe.    

clear all
close all
clc

%% Add multishape aligment functions to the search path 
disp('Load shape alignment library files.');
run ../lib/init_lib.m;

%% Load the fish data 
load ../data/cdfhc_data2D_fish.mat;
disp('Loaded fish data')
X = fish2d; % noiseless data
%% Add gaussian noise to the shapes
X_noise = X;
for iSh = 1:length(X) 
  X_noise{iSh} = X{iSh}(1:2:end, :) + 0.05*randn(size(X{iSh}(1:2:end, :)));
end
%% Center all shapes with respect to the origin
X = center_shapes(X);
X_noise = center_shapes(X_noise);

%% Compute TPS matrix of X vs X_noise
for ob = 1:length(X)
    K_XX_noise{ob} = compute_tps_K(X{ob}, X_noise{ob});
    K_XX{ob} = compute_tps_K(X{ob}, X{ob});
    [Q1{ob}, Q2{ob}] =  compute_QR_fact([X{ob} ones(size(X{ob}, 1), 1)]);
end


%% Display set of centered shapes before registration
hfig = figure();
overlay_shapes(X);
title('Shapes before registration')
xlim([-4 4])
ylim([-2 2])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow, pause(0.1)
close(hfig);


hfig = figure();
overlay_shapes(X_noise);
title('Shapes before registration (Noise added)')
xlim([-4 4])
ylim([-2 2])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow, pause(0.1)
close(hfig);

% normalized Information Potential
options.disp_progress = 1;
options.init_sigma = 0.2;
options.final_sigma = 0.02;
options.anealing = 0.98;
options.max_iter = 150;

%% lambda = 0.1
options.regularization = 0.1;
[tXnmip_1, A_1, W_1] = IP_Norm_Diff_FixedPoint(X, options);
[tXnmip_noise_1, A_noise_1, W_noise_1] = IP_Norm_Diff_FixedPoint(X_noise, options);
hfig = figure(); 
overlay_shapes(tXnmip_1);
title('Normalized IP difference')
xlim([-4 4])
ylim([-2 2])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow, pause(0.1)
close(hfig);
hfig = figure();
overlay_shapes(tXnmip_noise_1);
title('Normalized IP difference (Noisy shapes)')
xlim([-4 4])
ylim([-2 2])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow, pause(0.1)
close(hfig);
% Mapping original points using the estimated transformation
% from noisy input data
hfig = figure();
tXnmip_rec_1 = tXnmip_1;
for ob = 1 : length(X) 
    tXnmip_rec_1{ob} = [X{ob} ones(size(X{ob}, 1), 1)]*A_noise_1{ob} + Q2{ob}*Q2{ob}'*K_XX_noise{ob}*W_noise_1{ob};
end
overlay_shapes(tXnmip_rec_1);
title('Normalized IP difference reconstruction')
xlim([-4 4])
ylim([-2 2])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow, pause(0.1)
close(hfig);


%% lambda = 0.01
options.regularization = 0.01;
[tXnmip_01, A_01, W_01] = IP_Norm_Diff_FixedPoint(X, options);
[tXnmip_noise_01, A_noise_01, W_noise_01] = IP_Norm_Diff_FixedPoint(X_noise, options);
hfig = figure(); 
overlay_shapes(tXnmip_01);
title('Normalized IP difference')
xlim([-4 4])
ylim([-2 2])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow, pause(0.1)
close(hfig);
hfig = figure();
overlay_shapes(tXnmip_noise_01);
title('Normalized IP difference (Noisy shapes)')
xlim([-4 4])
ylim([-2 2])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow, pause(0.1)
close(hfig);
% Mapping original points using the estimated transformation
% from noisy input data
hfig = figure();
tXnmip_rec_01 = tXnmip_01;
for ob = 1 : length(X) 
    tXnmip_rec_01{ob} = [X{ob} ones(size(X{ob}, 1), 1)]*A_noise_01{ob} + Q2{ob}*Q2{ob}'*K_XX_noise{ob}*W_noise_01{ob};
end
overlay_shapes(tXnmip_rec_01);
title('Normalized IP difference reconstruction')
xlim([-4 4])
ylim([-2 2])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow, pause(0.1)
close(hfig);



%% lambda = 0.001
options.regularization = 0.001;
[tXnmip_001, A_001, W_001] = IP_Norm_Diff_FixedPoint(X, options);
[tXnmip_noise_001, A_noise_001, W_noise_001] = IP_Norm_Diff_FixedPoint(X_noise, options);
hfig = figure(); 
overlay_shapes(tXnmip_001);
title('Normalized IP difference')
xlim([-4 4])
ylim([-2 2])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow, pause(0.1)
close(hfig);
hfig = figure();
overlay_shapes(tXnmip_noise_001);
title('Normalized IP difference (Noisy shapes)')
xlim([-4 4])
ylim([-2 2])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow, pause(0.1)
close(hfig);
% Mapping original points using the estimated transformation
% from noisy input data
hfig = figure();
tXnmip_rec_001 = tXnmip_001;
for ob = 1 : length(X) 
    tXnmip_rec_001{ob} = [X{ob} ones(size(X{ob}, 1), 1)]*A_noise_001{ob} + Q2{ob}*Q2{ob}'*K_XX_noise{ob}*W_noise_001{ob};
end
overlay_shapes(tXnmip_rec_001);
title('Normalized IP difference reconstruction')
xlim([-4 4])
ylim([-2 2])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow, pause(0.1)
close(hfig);



%% lambda = 0.0001
options.regularization = 0.0001;
[tXnmip_0001, A_0001, W_0001] = IP_Norm_Diff_FixedPoint(X, options);
[tXnmip_noise_0001, A_noise_0001, W_noise_0001] = IP_Norm_Diff_FixedPoint(X_noise, options);
hfig = figure(); 
overlay_shapes(tXnmip_0001);
title('Normalized IP difference')
xlim([-4 4])
ylim([-2 2])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow, pause(0.1)
close(hfig);
hfig = figure();
overlay_shapes(tXnmip_noise_0001);
title('Normalized IP difference (Noisy shapes)')
xlim([-4 4])
ylim([-2 2])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow, pause(0.1)
close(hfig);
% Mapping original points using the estimated transformation
% from noisy input data
hfig = figure();
tXnmip_rec_0001 = tXnmip_0001;
for ob = 1 : length(X) 
    tXnmip_rec_0001{ob} = [X{ob} ones(size(X{ob}, 1), 1)]*A_noise_0001{ob} + Q2{ob}*Q2{ob}'*K_XX_noise{ob}*W_noise_0001{ob};
end
overlay_shapes(tXnmip_rec_0001);
title('Normalized IP difference reconstruction')
xlim([-4 4])
ylim([-2 2])
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow, pause(0.1)
close(hfig);




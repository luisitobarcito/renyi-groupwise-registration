function [tX, A, W] = IP_Norm_Diff_FixedPoint(Xin, varargin)
% [tX, A, W] = IP_Norm_Diff_FixedPoint(X, options)
% Multi-shape registration using Normalized Information Potential differences
% This is a multi-shape non-rigid alignment algorithm
% It employes Thin-plate splines or radial basis function expansion
% for the non-rigid transofrmation and affine transofrmation for the
% rigid component. The algoirthm is based on a fixed point update rule
% derived from the gradient fo the normalized information potential difference.
% The algorithm is described in:
% "Group-wise point-set Registration Based on Renyi's Second Order
% Entropy"
%--------------------------------------------------------------------------
% Inputs
% Xin: A cell array of length N containing the shapes. Each shape is
%    represented by (M x d) matrix where M is the number of points and
%    d is the dimensionality of the shape.
%
% options: This is an optional parameter that specify custom options
%          for the algorithms. "options" is a struct with the following fields.
% options.nonrigid: String, 'TPS' or 'GaussianRBF' (default 'TPS')
% options.regularization: Tradeoff parameter for the regularization
%                        (default 0.001)
% options.init_sigma: Kernel bandwidth for the information potential
%                     estimator (default 0.1)
% options.final_sigma: Smallest value of kernel bandwidth for anealing
%                      (default 0.01)
% options.anealing: Exponential decay factor for the anealing the
%                   kernel bandwidth (default 0.97)
% options.max_iter: Number of fixed point iterations (default 300)
% options.disp_progress : plots_progress (default false)
% -------------------------------------------------------------------------
% Outputs
% tX: A cell array of length N containing the registered shapes. Each
%     shape of the same size as in X.
% A: A cell array of length N. Each cell contains a matrix of size
%    ((d+1) x d)the affine component of the
%    trasformation.
% W: A cell array of length N. Each cell contains a matrix (M x d)
%    the coefficients for the nonlinear deformations.

% Default parameters
default.nonrigid = 'TPS';
default.regularization = 0.001;
default.init_sigma = 0.1;
default.final_sigma = 0.01;
default.anealing = 0.97;
default.max_iter = 300;
default.disp_progress = false;
default.create_video = false;
default.video_file = 'ip_norm_video.avi';


%% Verify options and set algorithm parameters
if ~isempty(varargin)
    options = varargin{1};
    options = verifyFields(options, default);
else
    options = default;
end

% Number of shapes
objNum = length(Xin);
dim = size(Xin{1},2);

lmbd2 = options.regularization;
sig = options.init_sigma;
anealing = options.anealing;

itrNum = options.max_iter;

if options.disp_progress
    figh = figure;
    colmark = ['r'; 'b'; 'g'; 'm'; 'c'; 'y'; 'k'];
end

if options.create_video
    vid_wrt = VideoWriter(options.video_file);
    open(vid_wrt);
end


%% Extend x's to homogenous coordinates {x=[x,1]}
% and compute QR factorization 

X = cell(objNum,1);
Q1 = cell(objNum, 1);
Q2 = cell(objNum, 1);

for ob = 1:objNum
    X{ob} = [Xin{ob}, ones(length(Xin{ob}),1)];
    [Q1{ob}, Q2{ob}] = compute_QR_fact(X{ob});
end
%% Identify lengths and dimensions of the objects
lenX = zeros(objNum,1);
for ob = 1:objNum
    lenX(ob) = size(X{ob},1);
end
totalX = sum(lenX);
%% Create 1 vectors of lenX
OneX = cell(objNum,1);
for i = 1:objNum
    OneX{i} = ones(lenX(i),1);
end

%% Compute K on X point set
K = cell(objNum,1);
for ob = 1:objNum
    if strcmp(options.nonrigid, 'TPS')
      K{ob} = compute_tps_K(X{ob}(:,1:dim), X{ob}(:,1:dim));
    elseif strcmp(options.nonrigid, 'Gaussian_RBF')
      K{ob} = compute_gauss_K(X{ob}, X{ob}, .1, false);
    else
      error('nonrigid basis is not supported');
    end
    % Reduce K to lie on the kernel of the span of X
    K{ob} = Q2{ob}*Q2{ob}'*K{ob}*Q2{ob}*Q2{ob}';
    %K{ob} = K{ob}*Q2{ob}*Q2{ob}';
end




%% Initial guess for transformation paramaters
% Intialize the transformation matrices A with extended row (shift)
A = cell(objNum,1);
for ob = 1:objNum
    A{ob} = [eye(dim); zeros(1,dim)];
end
% Intialize the transformation matrices W
W = cell(objNum,1);
for ob = 1:objNum
    W{ob} = Q2{ob}*Q2{ob}'*pinv(K{ob})*(0.0001*randn(lenX(ob), dim)); %???
end

tXA = cell(objNum,1);
tKW = cell(objNum,1);
tX = cell(objNum,1);
CX = cell(objNum,1);
Ctot = zeros(dim);
Xtot_mean = zeros(1, dim);
for ob = 1:objNum
    tXA{ob} = X{ob}*A{ob};
    tKW{ob} =  K{ob}*W{ob};
    tX{ob} = tXA{ob} + tKW{ob};
    CX{ob} = tX{ob}'*tX{ob}/lenX(ob) - (tX{ob}'*OneX{ob})*(OneX{ob}'*tX{ob})/(lenX(ob)^2);
%    Ctot = Ctot + CX{ob}*(lenX(ob)/totalX);
    Ctot = Ctot + tX{ob}'*tX{ob};
    Xtot_mean = Xtot_mean + sum(tX{ob},1);
end
Ctot = Ctot/totalX - (Xtot_mean'*Xtot_mean)/(totalX^2); 


if options.disp_progress
    figure(figh), hold off;
    for ob = 1:objNum
        colidx = mod(ob-1, length(colmark)) + 1;
        plot(tX{ob}(:,1), tX{ob}(:,2), ['o',colmark(colidx)], 'MarkerSize',3, 'MarkerFaceColor', colmark(colidx), 'LineWidth', 3); hold on;
    end
    set(gca, 'DataAspectRatio', [1 1 1]);
    drawnow;
    pause(0.5);
    xlimits = get(gca, 'xlim');
    ylimits = get(gca, 'ylim');
end

cipX = zeros(objNum);
G = cell(objNum, objNum);
trCX = zeros(objNum,1);
Jcost = zeros(itrNum+1,1);
for itr = 1:itrNum
    %% ########################################################
    % Compute Gaussians (Y,tX) and (tX,tX), where tX = X*A + K*W
    % and compute cross-informations potentials and Covariance
    % matrices determinants and inverses
    % will compute only upper triangle
    for ob1 = 1:objNum          %ob1 = 1:obN
        for ob2 = ob1:objNum    %ob2 = ob1:obN
            G{ob1,ob2} = compute_gauss_K(tX{ob1}, tX{ob2}, sig, true);
            cipX(ob1, ob2) = OneX{ob1}'*G{ob1,ob2}*OneX{ob2}/(lenX(ob1)*lenX(ob2));
            G{ob2,ob1} = G{ob1,ob2}';
            cipX(ob2, ob1) = cipX(ob1, ob2);
        end
        trCX(ob1) = trace(CX{ob1});
    end
    trCtot = trace(Ctot);
    ipXtot = ((lenX'/totalX)*cipX*(lenX/totalX));
    normIpXtot = ipXtot/sqrt(trCtot); % Normalized IP of the Union
    normIpAvg  = ((lenX.^2)/(lenX'*lenX))'*(diag(cipX)./sqrt(trCX));% Normalized average IP
    
    %% Compute objective value
    Jcost(itr+1) = normIpAvg - normIpXtot;
    
    
    
    %% ##### Fixed point update A and W for each object#####
    for ob = 1:objNum
        
        %% Information potentials Covariance for each shape
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Information Potential
        dIPcoeff = (2/((lenX'*lenX)*sqrt(trCX(ob))*(sig^2))); % just a coeff
        Gkk_A_X = dIPcoeff*X{ob}'*(G{ob, ob})*X{ob};
        Gkk_A_K = dIPcoeff*X{ob}'*(G{ob, ob})*K{ob};
        diagGkk_A_X = dIPcoeff*X{ob}'*(diag(G{ob, ob}*OneX{ob}))*X{ob};
        diagGkk_A_K = dIPcoeff*X{ob}'*(diag(G{ob, ob}*OneX{ob}))*K{ob};
        
        Gkk_W_X = dIPcoeff*K{ob}'*(G{ob, ob})*X{ob};
        Gkk_W_K = dIPcoeff*K{ob}'*(G{ob, ob})*K{ob};
        diagGkk_W_X = dIPcoeff*K{ob}'*(diag(G{ob, ob}*OneX{ob}))*X{ob};
        diagGkk_W_K = dIPcoeff*K{ob}'*(diag(G{ob, ob}*OneX{ob}))*K{ob};
        %%%%% Covariance
        %dCkcoeff = (cipX(ob, ob)*(lenX(ob)^2)/((lenX'*lenX)*sqrt(trCX(ob))));
	dCkcoeff = (cipX(ob, ob)*(lenX(ob)^2)/((lenX'*lenX)*trCX(ob)*sqrt(trCX(ob))));
        Ck_A_X = (dCkcoeff/lenX(ob))*X{ob}'*X{ob};
        Ck_A_K = (dCkcoeff/lenX(ob))*X{ob}'*K{ob};
        Onek_A_X = (dCkcoeff/(lenX(ob)^2))*((X{ob}'*OneX{ob})*(OneX{ob}'*X{ob}));
        Onek_A_K = (dCkcoeff/(lenX(ob)^2))*((X{ob}'*OneX{ob})*(OneX{ob}'*K{ob}));
        
        Ck_W_X = (dCkcoeff/lenX(ob))*K{ob}'*X{ob};
        Ck_W_K = (dCkcoeff/lenX(ob))*K{ob}'*K{ob};
        Onek_W_X = (dCkcoeff/(lenX(ob)^2))*((K{ob}'*OneX{ob})*(OneX{ob}'*X{ob}));
        Onek_W_K = (dCkcoeff/(lenX(ob)^2))*((K{ob}'*OneX{ob})*(OneX{ob}'*K{ob}));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %% Cross-information potentials and Covariance for the Union of Shapes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Cross-Information Potential
        Gkl_A_tXA = zeros(size(X{ob},2), size(tXA{ob}, 2));
        Gkl_A_tKW = zeros(size(X{ob},2), size(tKW{ob}, 2));
        diagGkl_A_X = zeros(size(X{ob},2), size(X{ob}, 2));
        diagGkl_A_K = zeros(size(X{ob},2), size(K{ob}, 2));
        
        Gkl_W_tXA = zeros(size(K{ob},2), size(tXA{ob}, 2));
        Gkl_W_tKW = zeros(size(K{ob},2), size(tKW{ob}, 2));
        diagGkl_W_X = zeros(size(K{ob},2), size(X{ob}, 2));
        diagGkl_W_K = zeros(size(K{ob},2), size(K{ob}, 2));
        %%%%% Covariance terms
        sumOneXtXA = zeros(1, size(tXA{ob}, 2));
        sumOneXtKW = zeros(1, size(tKW{ob}, 2));
        
        dCIPcoeff =  (2/((totalX^2)*(sig^2)*sqrt(trCtot)));
	for ob2 = 1:objNum
            Gkl_A_tXA = Gkl_A_tXA + dCIPcoeff*X{ob}'*(G{ob, ob2})*tXA{ob2};
            Gkl_A_tKW = Gkl_A_tKW + dCIPcoeff*X{ob}'*(G{ob, ob2})*tKW{ob2};
            diagGkl_A_X = diagGkl_A_X  + dCIPcoeff*X{ob}'*(diag(G{ob,ob2}*OneX{ob2}))*X{ob};
            diagGkl_A_K = diagGkl_A_K +  dCIPcoeff*X{ob}'*(diag(G{ob,ob2}*OneX{ob2}))*K{ob};
            
            Gkl_W_tXA = Gkl_W_tXA +  dCIPcoeff*K{ob}'*(G{ob, ob2})*tXA{ob2};
            Gkl_W_tKW = Gkl_W_tKW +  dCIPcoeff*K{ob}'*(G{ob, ob2})*tKW{ob2};
            diagGkl_W_X = diagGkl_W_X + dCIPcoeff*K{ob}'*(diag(G{ob,ob2}*OneX{ob2}))*X{ob};
            diagGkl_W_K = diagGkl_W_K +  dCIPcoeff*K{ob}'*(diag(G{ob,ob2}*OneX{ob2}))*K{ob};
            
            sumOneXtXA = sumOneXtXA + OneX{ob2}'*tXA{ob2};
            sumOneXtKW = sumOneXtKW + OneX{ob2}'*tKW{ob2};
        end
	%dCtotcoeff = (ipXtot/sqrt(trCtot));
        dCtotcoeff = (ipXtot/(trCtot*sqrt(trCtot)));
        
	Ctot_A_X = (dCtotcoeff/totalX)*X{ob}'*X{ob};
        Ctot_A_K = (dCtotcoeff/totalX)*X{ob}'*K{ob};
        Onekl_A_tXA = (dCtotcoeff/(totalX^2))*((X{ob}'*OneX{ob})*(sumOneXtXA));
        Onekl_A_tKW = (dCtotcoeff/(totalX^2))*((X{ob}'*OneX{ob})*(sumOneXtKW));
        
        Ctot_W_X = (dCtotcoeff/totalX)*K{ob}'*X{ob};
        Ctot_W_K = (dCtotcoeff/totalX)*K{ob}'*K{ob};
        Onekl_W_tXA = (dCtotcoeff/(totalX^2))*((K{ob}'*OneX{ob})*(sumOneXtXA));
        Onekl_W_tKW = (dCtotcoeff/(totalX^2))*((K{ob}'*OneX{ob})*(sumOneXtKW));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %% CALCULATE  new  A  &  W
        AnumIP = (Gkl_A_tXA + Gkl_A_tKW - (Gkk_A_K - diagGkk_A_K + diagGkl_A_K)*W{ob} - (Gkk_A_X - diagGkk_A_X)*A{ob});
        AnumCov = (Ctot_A_X*A{ob} + Ctot_A_K*W{ob} - Onekl_A_tXA - Onekl_A_tKW - (Ck_A_X - Onek_A_X)*A{ob} - (Ck_A_K - Onek_A_K)*W{ob});
        Aden =   diagGkl_A_X;
%	AnumCov = (Ctot_A_X*A{ob} + Ctot_A_K*W{ob} - Onekl_A_tXA - Onekl_A_tKW - (0 - Onek_A_X)*A{ob} - (Ck_A_K - Onek_A_K)*W{ob});
%         Aden =   diagGkl_A_X - Ck_A_X;
        
        WnumIP = (Gkl_W_tXA + Gkl_W_tKW - (Gkk_W_X - diagGkk_W_X + diagGkl_W_X)*A{ob} - (Gkk_W_K  - diagGkk_W_K)*W{ob});
        WnumCov = (Ctot_W_X*A{ob} + Ctot_W_K*W{ob} - Onekl_W_tXA - Onekl_W_tKW - (Ck_W_X - Onek_W_X)*A{ob} - (Ck_W_K - Onek_W_K)*W{ob});
        Wden =    diagGkl_W_K + 2*lmbd2*K{ob};
%         WnumCov = (Ctot_W_X*A{ob} +  Ctot_W_K*W{ob} - Onekl_W_tXA - Onekl_W_tKW - (Ck_W_X - Onek_W_X)*A{ob} - (0 - Onek_W_K )*W{ob});
%         Wden =    diagGkl_W_K - Ck_W_K  + lmbd2*K{ob};
        
        
        newA{ob} = pinv(Aden)*(AnumIP - AnumCov);
        newW{ob} = pinv(Wden)*(WnumIP - WnumCov);
        
        
    end
    %#####end of A and W updates for each object#####
    
    %% Update shapes and covariances 
    Ctot = zeros(dim);
    Xtot_mean = zeros(1, dim);
    for ob = 1:objNum
	A{ob} = newA{ob};
	W{ob} = Q2{ob}*Q2{ob}'*newW{ob};
%	W{ob} = newW{ob};
        tXA{ob} = X{ob}*A{ob};
        tKW{ob} =  K{ob}*W{ob};
        tX{ob} = tXA{ob} + tKW{ob};
	CX{ob} = tX{ob}'*tX{ob}/lenX(ob) - (tX{ob}'*OneX{ob})*(OneX{ob}'*tX{ob})/(lenX(ob)^2);
%	Ctot = Ctot + CX{ob}*(lenX(ob)/totalX);
	Ctot = Ctot + tX{ob}'*tX{ob};
	Xtot_mean = Xtot_mean + sum(tX{ob},1);
    end
    Ctot = Ctot/totalX - (Xtot_mean'*Xtot_mean)/(totalX^2); 

    
    
    if options.disp_progress
        if(~mod(itr, options.disp_progress))
%             figure(figh), 
            hold off;
            for ob = 1:objNum
                colidx = mod(ob-1, length(colmark)) + 1;
                plot(tX{ob}(:,1), tX{ob}(:,2), ['o',colmark(colidx)], 'MarkerSize',3, 'MarkerFaceColor', colmark(colidx), 'LineWidth', 3); hold on;
            end
            set(gca, 'DataAspectRatio', [1 1 1], 'xlim', xlimits, 'ylim', ylimits);
            title([num2str(itr), '    sigma ', num2str(sig)]);  
            pause(0.01);
            if options.create_video
                frame = getframe(gcf);
                writeVideo(vid_wrt,frame);
            end
        end
    end
    % disp(num2str(sum(sum(  compute_gauss_K(Y, X, 1) / (lenY*lenX)  ))));
    
    
    if(~mod(itr,1))
        sig = max(sig*anealing, options.final_sigma);
    end
    
end

if options.create_video
    close(vid_wrt);
end

for ob = 1:objNum
    % Return transformed shapes
    tX{ob} = X{ob}*A{ob} + K{ob}*W{ob};
    % Project W to the orthogonal space of X
    W{ob} = Q2{ob}*Q2{ob}'*W{ob};
end
if exist('figh')
    close(figh);
end
end
%--------------------------------------------------------------------------



function options = verifyFields(options, default)
% Performs parameter verification of the options struct

% Complete unset fields with defaults;
default_fields = fieldnames(default);
for iFld = 1 : length(default_fields)
    if ~isfield(options, default_fields{iFld})
        options.(default_fields{iFld}) = default.(default_fields{iFld});
    end
end

if ~any(strcmp(options.nonrigid, {'TPS', 'Gaussian_RBF'}))
    warning('Invalid basis function type, Setting to default');
    options.nonrigid = default.nonrigid;
end

if options.regularization < 0
    warning('Regularization must be a nonnegative scalar. Setting to default');
    options.regularization = default.regularization;
end

if options.init_sigma <= 0
    warning('Initial bandwidth must be greater than zero. Setting to default');
    options.init_sigma = default.init_sigma;
end

if options.final_sigma <= 0
    warning('Final bandwidth must be greater than zero. Setting to default');
    options.final_sigma = default.final_sigma;
end

if options.anealing > 1 || options.anealing <=0
    warning('Anealing rate must be in the range (0, 1]. Setting to default');
    options.anelaing = default.anealing;
end

if options.max_iter < 10
    warning('Maximum number of iterations is too small. Setting to default');
    options.max_iter = default.max_iter;
elseif options.max_iter > 10000
    warning('Maximun number of itearations cannot be larger than 10000. max_iter will be set to 10000');
    options.max_iter = 10000;
end

if options.disp_progress < 0
    warning('Invalid value for display intervals. Setting display to false');
    options.disp_progress = false;
elseif options.disp_progress > options.max_iter
    warning('Value for display intervals is larger than the maximun number of iterations. No progress will be displayed');
end

if options.disp_progress == false && options.create_video == true 
    options.create_video = false;
    warning('Display intervals were not set. No video will be created');
end
end

%%
%Information Potential (IP) UPDATE AW using TPS
%--------------------------------------------------------------------------

%######################
%TPS
%----------------------
%initialization
%lambda = 1;
%sigma  = 1;
%----------------------
%annealing rate
%lambda = lambda*0.93;
%sigma  = sigma*0.98;
%----------------------
%iterations = 300
%######################


%v2 - change update of X - do not update X, just the A and W (based on old values)


function [tX, A, W] = Holder_FixedPoint(objects, varargin)
% [tX, A, W] = Holder_FixedPoint(X, options)
% Multi-shape registration using Holder's inequality
% This is a multi-shape non-rigid alignment algorithm
% It employes Thin-plate splines or radial basis function expansion
% for the non-rigid transofrmation and affine transofrmation for the
% rigid component. The algoirthm is based on a fixed point update rule
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

%% Verify options and set algorithm parameters
if ~isempty(varargin)
    options = varargin{1};
    options = verifyFields(options, default);
else
    options = default;
end

    objNum = length(objects);
    lmbd = options.regularization;
    sig = options.init_sigma;
    anealing = options.anealing;
    
    itrNum = options.max_iter;
    if options.disp_progress
      figh = figure;
      colmark = ['r'; 'b'; 'g'; 'm'; 'c'; 'y'; 'k'; 'r'; 'b'; 'g'; 'm'; 'c'; 'y'];
    end

    %extend x's to homogenous coordinates {x=[x,1]}
    X = cell(objNum,1);
    Q1 = cell(objNum, 1);
    Q2 = cell(objNum, 1);

    for ob = 1:objNum
        X{ob} = [objects{ob}, ones(length(objects{ob}),1)];
	[Q1{ob}, Q2{ob}] = compute_QR_fact(X{ob});
    end

    %identify lengths and dimensions of the objects
    lenX = zeros(objNum,1);
    OneX = cell(objNum,1);
    for ob = 1:objNum
        lenX(ob) = size(X{ob},1);
        OneX{ob} = ones(lenX(ob),1);
    end
    dim = size(X{1},2);
    
    %compute K on X point set
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
     % K{ob} = K{ob}*Q2{ob}*Q2{ob}';
    end
        
    %intialize the transformation matrices A and W, and lambda
    A = cell(objNum,1);
    for ob = 1:objNum
        A{ob} = eye(dim);
    end
    newA = cell(objNum,1);
    
    %intialize the transformation matrices W
    W = cell(objNum,1);
    for ob = 1:objNum
        W{ob} = Q2{ob}*Q2{ob}'*pinv(K{ob})*(0.0001*randn(lenX(ob), dim)); %???
    end
    newW = cell(objNum,1);
    
    
    aI = 1e-5*eye(dim);
    wI = cell(objNum,1);
    for ob = 1:objNum
    %     wI{ob} = 1e-5*eye(lenX(ob));
        wI{ob} = eye(lenX(ob));
    end

    tX = cell(objNum,1);
    for ob = 1:objNum
        tX{ob} = X{ob}*A{ob};
    end
    
    if options.disp_progress
      figure(figh), hold off;
      for ob = 1:objNum
	plot(tX{ob}(:,1), tX{ob}(:,2), ['o',colmark(ob)], 'MarkerSize',10, 'MarkerFaceColor', colmark(ob), 'LineWidth', 3); hold on;
      end
      set(gca, 'DataAspectRatio', [1 1 1]);
      drawnow;
      pause(0.5);
    end
 
    for itr = 1:itrNum
        
      %compute Guassians (Y,tX) and (tX,tX), where tX = X*A + K*W
      G = cell(objNum, objNum);   %will compute only upper triangle
      for ob1 = 1:objNum          %ob1 = 1:obN
        for ob2 = 1:objNum      %ob2 = 1:obN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          G{ob1,ob2} = compute_gauss_K(tX{ob1}, tX{ob2}, sig);
          G{ob1,ob2} = G{ob1,ob2} / sum(G{ob1,ob2}(:));
        end
      end
      
      %#####compute A for each object#####
      for ob = 1:objNum
          
        %FIRST TERM
        %calculate term Gxyz*X'*X
        %calculate and sum up all terms Gxyz*X'*Y*Ay
        GXYZ_X_YAs = zeros(dim);
        GXYZ_X_KWs = zeros(dim);
        GXYZ_K_YAs = zeros(lenX(ob),dim);
        GXYZ_K_KWs = zeros(lenX(ob),dim);
        
        GXYZ_X_X = zeros(dim);
        GXYZ_X_K = zeros(dim,lenX(ob));
        GXYZ_K_X = zeros(lenX(ob),dim);
        GXYZ_K_K = zeros(lenX(ob));
        for ob2 = [1:ob-1, ob+1:objNum]
          GXYZ_X_YAs = GXYZ_X_YAs + X{ob}'*G{ob,ob2}*X{ob2}*A{ob2};
          GXYZ_X_KWs = GXYZ_X_KWs + X{ob}'*G{ob,ob2}*K{ob2}*W{ob2};
          GXYZ_K_YAs = GXYZ_K_YAs + K{ob}*G{ob,ob2}*X{ob2}*A{ob2};
          GXYZ_K_KWs = GXYZ_K_KWs + K{ob}*G{ob,ob2}*K{ob2}*W{ob2};
          
          GXYZ_X_X = GXYZ_X_X + X{ob}'*diag(G{ob,ob2}*OneX{ob2})*X{ob};
          GXYZ_X_K = GXYZ_X_K + X{ob}'*diag(G{ob,ob2}*OneX{ob2})*K{ob};
          GXYZ_K_X = GXYZ_K_X + K{ob}*diag(G{ob,ob2}*OneX{ob2})*X{ob};
          GXYZ_K_K = GXYZ_K_K + K{ob}*diag(G{ob,ob2}*OneX{ob2})*K{ob};
        end
        
        %SECOND TERM

        GXXX_XX_XX = zeros(dim);
        GXXX_XX_KK = zeros(dim,lenX(ob));
        GXXX_KK_XX = zeros(lenX(ob),dim);
        GXXX_KK_KK = zeros(lenX(ob));
        for ob2 = [1:ob-1, ob+1:objNum]
          GXXX_XX_XX = GXXX_XX_XX + 2*X{ob}'*( diag(G{ob,ob}*OneX{ob}) - G{ob,ob} )*X{ob};
          GXXX_XX_KK = GXXX_XX_KK + 2*X{ob}'*( diag(G{ob,ob}*OneX{ob}) - G{ob,ob} )*K{ob};
          GXXX_KK_XX = GXXX_KK_XX + 2*K{ob}*( diag(G{ob,ob}*OneX{ob}) - G{ob,ob} )*X{ob};
          GXXX_KK_KK = GXXX_KK_KK + 2*K{ob}*( diag(G{ob,ob}*OneX{ob}) - G{ob,ob} )*K{ob};
        end
        
        
        %CALCULATE   A
        %First Term elements multiply by (objNum) {the power}
        newA{ob} = pinv(2*GXYZ_X_X) * ( 2*(GXYZ_X_YAs + GXYZ_X_KWs - GXYZ_X_K*W{ob}) + (GXXX_XX_XX*A{ob} + GXXX_XX_KK*W{ob}) );
        try 
            newW{ob} = pinv(2*GXYZ_K_K - lmbd*K{ob}) * ( 2*(GXYZ_K_YAs + GXYZ_K_KWs - GXYZ_K_X*A{ob}) + (GXXX_KK_XX*A{ob} + GXXX_KK_KK*W{ob}) );
        catch
            newW{ob} = (2*GXYZ_K_K - lmbd*K{ob}) \ ( 2*(GXYZ_K_YAs + GXYZ_K_KWs - GXYZ_K_X*A{ob}) + (GXXX_KK_XX*A{ob} + GXXX_KK_KK*W{ob}) );
        end
      end
      %#####end of A for each object#####
      
      

      for ob = 1:objNum
        if(mod(itr,2))
          A{ob} = newA{ob};
        else
          W{ob} = Q2{ob}*Q2{ob}'*newW{ob};
        end
        tX{ob} = X{ob}*A{ob} + K{ob}*W{ob};
      end
      
      % %         if(itr<300), X = X*A;
      % %         else        X = bsxfun(@plus, X*A, K*W);
      % %         end
      if options.disp_progress
        if(~mod(itr, options.disp_progress))
          figure(figh), hold off;
          for ob = 1:objNum
            plot(tX{ob}(:,1), tX{ob}(:,2), ['o',colmark(ob)], 'MarkerSize',10, 'MarkerFaceColor', colmark(ob), 'LineWidth', 3); hold on;
          end
          set(gca, 'DataAspectRatio', [1 1 1]);
          title([num2str(itr), '    sigma ', num2str(sig)]);
        end
        pause(0.01);
      end
      

      
      % disp(num2str(sum(sum(  compute_gauss_K(Y, X, 1) / (lenY*lenX)  ))));
      
      
      if(~mod(itr,1))
	%           lmbd = lmbd*0.95;
        sig = max(sig*anealing, options.final_sigma);
      %             disp(['lambda   ', num2str(lmbd), '      sigma     ', num2str(sig)]);
      end
    % pause(0.5);
    end

    %    disp(['lambda ', num2str(lmbd), '    sigma ', num2str(sig)]);
    
    
    for ob = 1:objNum
      A{ob} = newA{ob};
      W{ob} = Q2{ob}*Q2{ob}'*newW{ob};
      tX{ob} = X{ob}*A{ob} + K{ob}*W{ob};
    end
    if exist('figh')
      close(figh);
    end

    
end
%--------------------------------------------------------------------------


%%

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

end






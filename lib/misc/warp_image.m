function im_w = warp_image(im, params)

% Default parameters
default.nonrigid = 'TPS';
default.size = size(im);

%% Verify given parameters
options = verifyFields(params, default);
A = options.A;
W = options.W;

% set of centers for the nonrigid tranformation
Xc = params.centers;
% create input image lattice
xg = 1:size(im, 2);
yg = 1:size(im, 1);
Xg = [reshape(repmat(xg(:)',length(yg), 1), [],1), reshape(repmat(yg(:),1, length(xg)), [],1)];
Xg = bsxfun(@minus, Xg, params.origin)*params.scale;
%% Extend x's to homogenous coordinates {x=[x,1]}
% and compute QR factorization
Xc = [Xc, ones(length(Xc),1)];
Xg = [Xg, ones(length(Xg),1)];
[~, Q2] = compute_QR_fact(Xc);
%% Compute K on X point set
if strcmp(options.nonrigid, 'TPS')
    Kcc = compute_tps_K(Xc(:,1:2), Xc(:,1:2));
    Kgc = compute_tps_K(Xg(:,1:2), Xc(:,1:2));
elseif strcmp(options.nonrigid, 'Gaussian_RBF')
    Kcc = compute_gauss_K(Xc, Xc, .1, false);
    Kgc = compute_tps_K(Xg, Xc, .1, false);
    
else
    error('nonrigid basis is not supported');
end
% Reduce K to lie on the kernel of the span of X
Kcc_pinv = pinv(Kcc);
Kcc =  (Q2*Q2')*Kcc*(Q2*Q2');


tXA = Xg*A;
tKW =  Kgc*Kcc_pinv*Kcc*W;
tXg = tXA + tKW;

tXcA = Xc*A;
tKccW =  Kcc*W;
tXcg = tXcA + tKccW;

%% output image grid

x_min = options.x_min;
x_max = options.x_max;

y_min = options.y_min;
y_max = options.y_max;



txg = linspace(x_min, x_max, options.size(2));
tyg = linspace(y_min, y_max, options.size(1));
Yg = [reshape(repmat(txg(:)', length(tyg), 1), [],1), reshape(repmat(tyg(:), 1, length(txg)), [],1)];

idx = knnsearch(tXg, Yg);
im = im(end:-1:1, :,:);

im_w = zeros(options.size);
for iCh = 1:size(im, 3)
    im_w_ch = im_w(:,:,iCh);
    im_ch = im(:,:, iCh);
    im_w_ch(:) = im_ch(idx);
    im_w(:,:,iCh) = im_w_ch;
end
im_w = im_w(end:-1:1,:,:);
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

if ~any(strcmp(options.nonrigid, {'TPS', 'GaussianRBF'}))
    warning('Invalid basis function type, Setting to default');
    options.nonrigid = default.nonrigid;
end

assert(isfield(options, 'centers'), 'centers of nonrigid tranformation must be supplied');
assert(isfield(options, 'W'), 'nonrigid tranformation matrix W must be supplied');
assert(isfield(options, 'A'), 'Affine tranformation matrix A must be supplied');


end



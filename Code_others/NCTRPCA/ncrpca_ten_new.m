function [D, I] = ncrpca_ten_new(T, true_r, params)
% function [D, I] = ncrpca_ten_new(T, true_r, EPS,MAX_ITER,EPS_S,incoh,TOL)

T = tensor(T);

if nargin < 3
    params = struct();
end

if ~isfield(params, 'TOL')
    TOL = 1e-1;
else
    TOL = params.TOL;
end

if ~isfield(params, 'incoh')
    incoh = 1;
else
    incoh = params.incoh;
end

if ~isfield(params, 'EPS_S')
    EPS_S = 1e-3;
else
    EPS_S = params.EPS_S;
end

if ~isfield(params, 'MAX_ITER')
    MAX_ITER = 100;
else
    MAX_ITER = params.MAX_ITER;
end

if ~isfield(params, 'EPS')
    EPS = 1e-3;
else
    EPS = params.EPS;
end

% if nargin < 7, TOL = 1e-1; end
% if nargin < 6, incoh = 1; end
% if nargin < 5, EPS_S = 1e-3; end
% if nargin < 4, MAX_ITER = 51; end
% if nargin < 3, EPS = 1e-3; end

frob_err(1) = inf;
n = size(T);
n1 = n(1);
n2 = n(2);
n3 = n(3);
t = 1;
idx = [];

if ~isfield(params, 'thresh_const')
    thresh_const = 1e3; % threshold constant: can be tuned depending on incoherence
else
    thresh_const = params.thresh_const;
end
if ~isfield(params, 'thresh_red')
    thresh_red = 0.9; % parameter to reduce the threshold constant: can be tuned
else
    thresh_red = params.thresh_red;
end

r_hat = 1; % initial rank for stagewise algorithm
L_t = tensor(zeros(size(T)));
Sig_t = cp_als(T,1,'printitn',0);
Sig_t = Sig_t.lambda(1);
D_t = T-L_t;
thresh = thresh_const*Sig_t/sqrt(n1*n2*n3);
idx = unique([find(abs(double(D_t)) > thresh); idx]);
S_t = tensor(zeros(size(T)));
S_t(idx) = D_t(idx); % initial thresholding
if max(idx(:))==0
    idx = [];
end
while frob_err(t)/norm(T)>=EPS && t<MAX_ITER % convergence check
    if ~mod(t, 10) % check progress
        fprintf('Iter no. %d\n', t);
    end
    t = t+1;
    L_t = cp_als(T-S_t, r_hat, 'printitn', 0);
    %sig_min = L_t.lambda(r_hat);
    L_t = tensor(L_t);
    D_t = double(T-L_t);
    thresh = (thresh_const/sqrt(n1*n2*n3));%*sig_min;%*w(r_hat, r_hat);
    idx = unique([find(abs(D_t) > thresh); idx]);
    S_t(idx) = D_t(idx);
    frob_err(t) = norm(T-(L_t+S_t));
    if ((frob_err(t-1)-frob_err(t))/frob_err(t-1) <= TOL) && r_hat<true_r
        r_hat = r_hat+1; % use this for incrementally updating rank by 1
%         sig_t = lansvd(M-S_t, true_r, 'L'); % svd function from propack
%         ratio_sig = sig_t(r_hat+1:end)./[sig_t(r_hat+2:end); sig_t(end)];
%         [~, mx_idx] = max(ratio_sig);
%         r_hat = r_hat+mx_idx; % update rank for the next stage
    elseif ((frob_err(t-1)-frob_err(t))/frob_err(t-1) <= TOL) && r_hat==true_r
        thresh_const = thresh_const*thresh_red; % tune threshold
    end
%     if ((frob_err(t-1)-frob_err(t))/frob_err(t-1) <= TOL) && r_hat<=1.5*min([n1,n2,n3])
%         r_hat = r_hat+1;
%     elseif((frob_err(t-1)-frob_err(t))/frob_err(t-1) <= TOL) && r_hat==1.5*min([n1,n2,n3])%true_r
%         thresh_fact = thresh_fact*thresh_mult;
%     end
end
S_t = double(S_t);
S_t(abs(S_t)<EPS_S) = 0;
S_t = tensor(S_t);
iters = length(frob_err)-1;

D.E = double(S_t);
D.L = double(L_t);
I.iter = iters;
I.err = frob_err(2:end);

end

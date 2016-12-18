function [T_t,L_t,S_t,iters,frob_err] = ls_alt_min_ten_ktensor(T,true_r,r_hat,EPS,MAX_ITER,EPS_S)
t = 1;
frob_err(1) = inf;
T_t = T;
n = size(T);
n1 = n(1);
n2 = n(2);
n3 = n(3);
idx = [];
thresh_fact = 1000;
thresh_mult = .98;
TOL = 1e-4;
%r_hat = 1;
while frob_err(t)/norm(T)>=EPS && t<MAX_ITER % convergence check
    if ~mod(t, 10) % check progress
        fprintf('Iter no. %d\n', t);
    end
    t = t+1;
    L_t = cp_als(T_t, r_hat, 'printitn', 0);
    %sig_min = L_t.lambda(r_hat);
    %L_t = tensor(L_t);
    D_t = T_t-L_t;
    thresh = (thresh_fact/sqrt(n1*n2*n3))^(1/3);%*sig_min;%*w(r_hat, r_hat);
    idx = unique([find(abs(D_t.U{1}) > thresh); idx]);
    T_t.U{1}(idx) = L_t.U{1}(idx);
    idx = unique([find(abs(D_t.U{2}) > thresh); idx]);
    T_t.U{2}(idx) = L_t.U{2}(idx);
    T_t.U{3} = [diag(ones(true_r,1)); zeros(n3-true_r,true_r)];
    S_t = T-T_t;
    frob_err(t) = norm(T.U{1}-(L_t.U{1}+S_t.U{1}),'fro');
    if ((frob_err(t-1)-frob_err(t))/frob_err(t-1) <= TOL)
        thresh_fact = thresh_fact*thresh_mult;
    end
%     if ((frob_err(t-1)-frob_err(t))/frob_err(t-1) <= TOL) && r_hat<=1.5*min([n1,n2,n3])
%         r_hat = r_hat+1;
%     elseif((frob_err(t-1)-frob_err(t))/frob_err(t-1) <= TOL) && r_hat==1.5*min([n1,n2,n3])%true_r
%         thresh_fact = thresh_fact*thresh_mult;
%     end
end
%S_t = double(S_t);
S_t.U{1}(abs(S_t.U{1})<EPS_S) = 0;
S_t.U{2}(abs(S_t.U{2})<EPS_S) = 0;
S_t.U{3}(abs(S_t.U{3})<EPS_S) = 0;
%S_t = tensor(S_t);
iters = length(frob_err)-1;
end

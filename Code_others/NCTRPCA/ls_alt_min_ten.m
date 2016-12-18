function [T_t,L_t,S_t,iters,frob_err] = ls_alt_min_ten(T, r_hat, EPS, MAX_ITER, EPS_S)
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
    L_t = tensor(L_t);
    D_t = double(T_t-L_t);
    thresh = (thresh_fact/sqrt(n1*n2*n3));%*sig_min;%*w(r_hat, r_hat);
    idx = unique([find(abs(D_t) > thresh); idx]);
    T_t(idx) = L_t(idx);
    S_t = T-T_t;
    frob_err(t) = norm(T-(L_t+S_t));
    if ((frob_err(t-1)-frob_err(t))/frob_err(t-1) <= TOL)
        thresh_fact = thresh_fact*thresh_mult;
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
end

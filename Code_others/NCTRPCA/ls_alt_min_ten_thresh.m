function [T_t, L_t, S_t, iters, frob_err] = ls_alt_min_ten_thresh(T, r_hat, EPS, MAX_ITER, EPS_S)
t = 1;
frob_err(1) = inf;
T_t = T; % initialize M_t as M, not as zeros(m, n);
idx = [];
thresh = max(max(max(double(T))));
thresh_rate = 0.6;
while frob_err(t)/norm(T)>=EPS && t<MAX_ITER % convergence check
    if ~mod(t, 10) % check progress
        fprintf('Iter no. %d\n', t);
    end
    if t>3 && frob_err(t)-frob_err(t-1)<1e-3
        thresh = thresh*thresh_rate;
    end
    t = t+1;
    L_t = tensor(cp_als(T_t, r_hat, 'printitn', 0));
    D_t = double(T_t-L_t);
    idx = find(abs(D_t) > thresh);
    T_t(idx) = L_t(idx);
    S_t = double(T-tensor(T_t));
    S_t(abs(S_t)<EPS_S) = 0;
    S_t = tensor(S_t);
    frob_err(t) = norm(T-(L_t+S_t));
end
iters = length(frob_err)-1; % no. of iters. of am = length(frob_err)-1
end

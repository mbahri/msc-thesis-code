function [ Uc, Ur, T, E, PlotsUc, PlotsUr ] = vbbrpca_kronfac( Y, T_in, E_in, Uc_in, Ur_in, r, MAXITER, O, E_ori)
%VBBRPCA_KRONFAC VBBRPCA with Kronecker factorization
addpath Bayes_tools/
addpath Tools/
addpath mmx_package/

% Crap for plotting
PlotsUc = cell(1, MAXITER);
PlotsUr = cell(1, MAXITER);

% MAXITER = 1;

% Constants about the data
[n, m] = size(Y(:,:,1));
Nobs = size(Y, 3);
% M = mean(T_in, 3);
M = zeros(r);

% Initialization of gamma, alpha, beta
a = 1e-5;
b = 1e-5;

Y2sum = sum(Y(:).^2);
scale2 = Y2sum / (Nobs*m*n);
scale = sqrt(scale2);

% gamma = gamrnd(a, b, 1, r);
gamma = scale*ones(1,r);

% gamma = scale*ones(r,1);        
beta = 1./scale2;
alpha = ones(n,m,Nobs)*scale;

% Work variables
T = T_in;
E = E_in;
Uc = Uc_in;
Ur = Ur_in;

% Temp vars
% mu_T = zeros(r^2, Nobs);

Sigma_Uc = scale*eye(r,r);
Sigma_Ur = scale*eye(r,r);

% We use the expectations of Uc'Uc and Ur'Ur a lot so we precompute
% them
%tic;
UcUc = Uc'*Uc + n*Sigma_Uc;
vUcUc = to_vec(UcUc);
%toc
%tic;
UrUr = Ur'*Ur + m*Sigma_Ur;
vUrUr = to_vec(UrUr);
%toc
    
% Main loop
for it = 1:MAXITER,
    
    tic;
    
    % Update covariance of E
    Sigma_E = 1 ./ (beta + alpha);
    
    % Update of En and Tn
    beta_alpha_ratio = ( beta * Sigma_E );
%     tic;
    Sigma_T = (eye(r^2) + beta*kron(UrUr, UcUc))^(-1);
%     t = toc;
%     ppup('Sigma_T', t);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     tic;
    [BT,CT] = sparse_krondecomp( Sigma_T, r, r, r, r );
%     t = toc;
%     fprintf('Factorized Sigma_T in %fs - relative error: %f\n', t, norm(Sigma_T - kron(BT,CT), 'fro') / norm(Sigma_T, 'fro') );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Update E first
%     tic;
    Temp = left_right_matrix_tensor_prod(Uc, Ur, T, 'nn', 'nt');
    E = beta_alpha_ratio .* (Y - Temp);
%     t = toc;
%     ppup('E', t);
    
    % Precompute Y - E, used multiple times after
    YmE = Y - E;
    
    % Update T
%     tic;
    Temp = left_right_matrix_tensor_prod(Uc, Ur, YmE, 'tn', 'nn');
    Temp = Temp + repmat(M, 1, 1, Nobs);
    for k=1:Nobs
        T(:,:,k) = from_vec( Sigma_T * to_vec(Temp(:,:,k)), r, r );
    end
%     t = toc;
%     ppup('T', t);
    
    %%%%%%%%%%% Update Uc and Ur
    gamma_diag = diag(gamma);
    % Precompute the expectations of T kron T - faster but NOT SUITABLE FOR
    % RANKS > 50
%     kronTpV = zeros(r^2, r^2, Nobs);
%     tic;
    V = get_V(Sigma_T, r);
%     [B,C] = sparse_krondecomp( V, r, r, r, r );
%     t = toc;
%     fprintf('Factorized V in %fs - relative error: %f\n', t, norm(V - kron(B,C), 'fro') / norm(V, 'fro') );
%     fprintf('Computed V in %fs\n', t);
%     for k=1:Nobs
%         kronTpV(:,:,k) = kron(T(:,:,k), T(:,:,k)) + V;
%     end
    % Precompute the difference between Y and E
    
    % Uc
%     tic;
    SumW = zeros(r);
%     SumR = zeros(r, n);
    Temp = mmx('mult', Ur, YmE, 'tt');
    SumR = sum(mmx('mult', T, Temp), 3);
    tUrCB = trace(UrUr*CT)*BT;
    for k=1:Nobs
        Tk = T(:,:,k);
%         SumW = SumW + from_vec( (kron(T(:,:,k), T(:,:,k)) + V) * vUrUr , r, r );
%         SumW = SumW + Tk*UrUr*Tk' + C*UrUr*B';
        SumW = SumW + Tk*UrUr*Tk' + tUrCB;
%         SumR = SumR + T(:,:,k) * Ur' * YmE(:,:,k)';
    end
    Sigma_Uc = beta*(gamma_diag + beta * SumW)^(-1);
    Uc = ( Sigma_Uc * SumR )';
%     t = toc;
%     ppup('Uc and Sigma_Uc', t);
    
    UcUc = Uc'*Uc + n*Sigma_Uc;
    vUcUc = to_vec(UcUc);
    
    % Ur
%     tic;
    SumW = zeros(r);
%     SumR = zeros(r, m);
    Temp = mmx('mult', Uc, YmE, 'tn');
    SumR = sum(mmx('mult', T, Temp), 3);
    tUcCB = trace(UcUc*CT)*BT;
    for k=1:Nobs
        Tk = T(:,:,k);
%         SumW = SumW + from_vec( (kron(T(:,:,k), T(:,:,k)) + V)' * vUcUc , r ,r );
%         SumW = SumW + Tk'*UcUc*Tk + C'*UcUc*B;
        SumW = SumW + Tk'*UcUc*Tk + tUcCB;
%         SumR = SumR + T(:,:,k)' * Uc' * YmE(:,:,k);
    end
    Sigma_Ur = beta*(gamma_diag + beta * SumW)^(-1);
    Ur = ( Sigma_Ur * SumR )';
%     t = toc;
%     ppup('Ur and Sigma_Ur', t);

    UrUr = Ur'*Ur + m*Sigma_Ur;
    vUrUr = to_vec(UrUr);
    
    % Update alpha
%     alpha = Nobs ./ (Nobs * Sigma_E + sum(E, 3).^2);
    alpha = 1 ./ (Sigma_E + E.^2);  % CHANGED FOR TEST MULTIPLE ALPHA
    
    % Update beta
%     tic;
    sum_of_norms = 0;
%     for k=1:Nobs
%         sum_of_norms = sum_of_norms + ...
%             expectation_norm_fac(Y(:,:,k), Uc, T(:,:,k), Ur, E(:,:,k), ...
%             UcUc, UrUr, B, C, Sigma_E(:,:,k)); % CHANGED FOR TEST MULTIPLE ALPHA
%     end
    for k=1:Nobs
        sum_of_norms = sum_of_norms + ...
            expectation_norm(Y(:,:,k), Uc, T(:,:,k), Ur, E(:,:,k), ...
            vUcUc, vUrUr, Sigma_T, V, Sigma_E(:,:,k)); % CHANGED FOR TEST MULTIPLE ALPHA
    end
%     t = toc;
%     ppup('beta', t);
    
    beta = Nobs*m*n / sum_of_norms;
    
    % Update gamma
%     gamma_b = 2*b ./ (2 + b * diag(UcUc + UrUr)); 
%     gamma = (a + (m + n)/2) ./ gamma_b;
%     tic;
    gamma = (2*a + m + n) ./ (2*b + diag(UcUc + UrUr));
%     t = toc;
%     ppup('gammas', t);
    
    t = toc;
    
    % Output useful information
    fprintf('[%d] | Max_err_e: %f - Max_err_o: %f - Max_err_tot: %f | beta: %f\n', ...
        it, ...
        max_rec_e(E, E_ori), ...
        max_rec_o(Uc,T,Ur,O), ...
        max_rec_f(Uc,T,Ur,E,Y), ...
        beta ...
    );
    fprintf('       Estimated rank | Uc: %d - Ur: %d\n', ...
        estim_rank(Uc), ...
        estim_rank(Ur)...
    );
    fprintf('       Sparsity | Max alpha: %f - Max gamma: %f | Time: %f\n', ...
        max(alpha(:)), ...
        max(gamma(:)), ...
        t ...
    );

    % stemplot(Uc,1,1);
    % waitforbuttonpress;
    % stemplot(Ur,1,1);
    % waitforbuttonpress;
    % meshplot(Uc,1,1);
    % waitforbuttonpress;
    % meshplot(Ur,1,1);
    % waitforbuttonpress;
    % meshplot(T(:,:,1),1,1);
    % waitforbuttonpress;
    % meshplot(E(:,:,1),1,1);
    % waitforbuttonpress;
    % subplot(1,3,1), imagesc(E_ori(:,:,1))
    % colorbar;
    % subplot(1,3,2), imagesc(E_in(:,:,1))
    % colorbar;
    % subplot(1,3,3), imagesc(E(:,:,1))
    % colorbar;
    % waitforbuttonpress;

    % fprintf('Continuing\n');
end

end

function [ M ] = get_traces(A)
    r = size(A, 3);
    M = zeros(r);
    
    for i=1:r
        M(i,i) = trace(A(:,:,i));
    end
end

function ppup(Q, t)
%     fprintf('Updated %s in %fs\n', Q, t);
end

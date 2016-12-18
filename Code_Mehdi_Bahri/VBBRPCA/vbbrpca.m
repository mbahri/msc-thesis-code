function [ Data, Info ] = vbbrpca( Y, r, MAXITER, Ori, Init)
%VBBRPCA VBBRPCA with alternative computation of Uc and Ur

% addpath Bayes_tools/
% addpath Tools/

% T_in, E_in, Uc_in, Ur_in, 

if nargin < 4 || ~isfield(Ori, 'O')
    O = ones(size(Y));
else
    O = Ori.O;
end

if nargin < 4 || ~isfield(Ori, 'E')
    E_ori = ones(size(Y));
else
    E_ori = Ori.E;
end

% Crap for plotting
PlotsUc = cell(1, MAXITER);
PlotsUr = cell(1, MAXITER);

% MAXITER = 1;

% Constants about the data
[n, m] = size(Y(:,:,1));
Nobs = size(Y, 3);

%% Initialization of gamma, alpha, beta
a_gamma0 = 1e-7;
b_gamma0 = 1e-5;
a_beta0  = 1e-6;
b_beta0  = 1e-6;
a_alpha0 = 1e-6;
b_alpha0 = 1e-6;
initVar = 1;

%% BABACAN STYLE INITIALIZATION
% Y2sum = sum(Y(:).^2);
% scale2 = Y2sum / (Nobs*m*n);
% scale = sqrt(scale2);
% 
% gamma = scale*ones(1,r);
% 
% beta = 1./scale2;
% alpha = ones(n,m,Nobs)*scale;

gamma = (a_gamma0+eps)/(b_gamma0+eps)*ones(r,1);
beta = (a_beta0+eps)/(b_beta0+eps);
alpha = initVar.^(-1)*ones(n,m,Nobs).*((a_alpha0+eps)/(b_alpha0+eps));

% Work variables that are given the intial values we get from the
% deterministic algorithm or from ML initalization, or random init.

if nargin < 5
    [D, ~] = rpca2d_l1(Y, 'maxiter', 0, 'r', r, 'time', 2);
    T = D.T;
    Uc = D.Uc;
    Ur = D.Ur;
    E = Y - wm_make_L( Uc, Ur, T );
else
    T = Init.T;
    Uc = Init.Uc;
    Ur = Init.Ur;
    E = Init.E;
end

% Sigma_Uc = scale*eye(r);
% Sigma_Ur = scale*eye(r);

Sigma_Uc = diag(gamma);
Sigma_Ur = diag(gamma);

% %% RANDOM INITIALIZATION ACCORDING TO THE PRIOR
% gamma = gamrnd(a, b, 1, r);
% 
% % If p(x) is prop to 1/x, then the CFD is log(x) so we can simulate a
% % sample by inverting the CFD
% beta = exp(rand());
% alpha = sqrt(exp(rand(n,m,Nobs)));
% 
% T = randn(r,r,Nobs);
% Uc = repmat(gamma, n, 1);
% Ur = repmat(gamma, m, 1);
% E = randn(n,m,Nobs) ./ alpha;
% 
% Sigma_Uc = diag(gamma);
% Sigma_Ur = diag(gamma);

%% Mean and initialization of the inner products
% We use the expectations of Uc'Uc and Ur'Ur a lot so we precompute
% them. NOTE: THIS IS PROBABLY WHERE THE PROBLEM IS
%tic;
UcUc = Uc'*Uc + n*Sigma_Uc;
vUcUc = to_vec(UcUc);
%toc
%tic;
UrUr = Ur'*Ur + m*Sigma_Ur;
vUrUr = to_vec(UrUr);
%toc

% M = mean(T, 3);
M = zeros(r);
    
%% Main loop
for it = 1:MAXITER
    
    
    % Update covariance of E
    Sigma_E = 1 ./ (beta + alpha);
    
    % Update of En and Tn
    beta_alpha_ratio = ( beta * Sigma_E );
    Sigma_T_i = eye(r^2) + beta*kron(UrUr, UcUc) ;
%     print_condition_number(Sigma_T_i);

%     tic;
    Sigma_T = Sigma_T_i^(-1);
%     toc
    
    % Tensorized operation
%     tic
%     E = beta_alpha_ratio .* (Y - ttm(tensor(T), {Uc, Ur}, [1 2]));
%     YmE = tensor(Y - E);
%     muT = repmat(M, 1, 1, Nobs) + beta*ttm(YmE, {Uc', Ur'}, [1 2]);
%     T = tensor(double(Sigma_T * tenmat(muT, 3)'), [r r Nobs]);
%     toc

%     tic
    for k=1:Nobs
        % Every E(i,j,k) can be computed at once by this vectorized
        % expression
        E(:,:,k) = beta_alpha_ratio(:,:,k) .* ( Y(:,:,k) - Uc*T(:,:,k)*Ur' );   % CHANGED FOR TEST MULTIPLE ALPHA
        
        % Compute mu1'
        T(:,:,k) = dlyap(-beta*(Uc'*Uc), Ur'*Ur, M + beta*Uc'*(Y(:,:,k) - E(:,:,k))*Ur);
    end
%     toc
    
    %%%%%%%%%%% Update Uc and Ur
    gamma_diag = diag(gamma);
    % Precompute the expectations of T kron T - faster but NOT SUITABLE FOR
    % RANKS > 50
%     kronTpV = zeros(r^2, r^2, Nobs);
    V = get_V(Sigma_T, r);
%     for k=1:Nobs
%         kronTpV(:,:,k) = kron(T(:,:,k), T(:,:,k)) + V;
%     end

    % Precompute the difference between Y and E
    YmE = Y - E;
    
    % Uc
    SumW = zeros(r);
    SumR = zeros(r, n);
    for k=1:Nobs
        SumW = SumW + from_vec( (kron(T(:,:,k), T(:,:,k)) + V) * vUrUr , r, r );
        SumR = SumR + T(:,:,k) * Ur' * YmE(:,:,k)';
    end
    Sigma_Uc_i = beta*(gamma_diag + beta * SumW);
%     print_condition_number(Sigma_Uc_i);
    Sigma_Uc = (Sigma_Uc_i)^(-1) ;
    Uc = (Sigma_Uc * SumR )';
    
    UcUc = Uc'*Uc + n.*Sigma_Uc;
    vUcUc = to_vec(UcUc);
    
    % Ur
    SumW = zeros(r);
    SumR = zeros(r, m);
    for k=1:Nobs
        SumW = SumW + from_vec( (kron(T(:,:,k), T(:,:,k)) + V)' * vUcUc , r ,r );
        SumR = SumR + T(:,:,k)' * Uc' * YmE(:,:,k);
    end
    Sigma_Ur_i = beta*(gamma_diag + beta * SumW);
%     print_condition_number(Sigma_Ur_i);
    Sigma_Ur = (Sigma_Ur_i)^(-1);
    Ur = (Sigma_Ur * SumR )';

    UrUr = Ur'*Ur + m.*Sigma_Ur;
    vUrUr = to_vec(UrUr);
    
    % Update alpha
%     alpha = Nobs ./ (Nobs * Sigma_E + sum(E, 3).^2);
    alpha = 1 ./ (Sigma_E + E.^2);  % CHANGED FOR TEST MULTIPLE ALPHA
    
    % Update beta
    sum_of_norms = 0;
    for k=1:Nobs
        sum_of_norms = sum_of_norms + ...
            expectation_norm(Y(:,:,k), Uc, T(:,:,k), Ur, E(:,:,k), ...
            UcUc, UrUr, Sigma_T, V, Sigma_E(:,:,k)); % CHANGED FOR TEST MULTIPLE ALPHA
    end
    
    beta = Nobs*m*n / sum_of_norms;
    
    % Update gamma
%     gamma_b = 2*b_gamma0 ./ (2 + b_gamma0 * diag(UcUc + UrUr)); 
%     gamma = (a_gamma0 + (m + n)/2) ./ gamma_b;
    gamma = (2*a_gamma0 + m + n) ./ (2*b_gamma0 + diag(UcUc + UrUr));

    t = toc;
    
    % Output useful information
    fprintf('[%d] | Max_err_e: %f - Max_err_o: %f - Max_err_tot: %f | beta: %f\n\n', ...
        it, ...
        max_rec_e(E, E_ori), ...
        max_rec_o(Uc,T,Ur,O), ...
        max_rec_f(Uc,T,Ur,E,Y), ...
        beta ...
    );
    fprintf('       Estimated rank | Uc: %d - Ur: %d\n\n', ...
        estim_rank(Uc), ...
        estim_rank(Ur)...
    );
    fprintf('       Sparsity | Max alpha: %f - Max gamma: %f | Time: %f\n\n', ...
        max(alpha(:)), ...
        max(gamma(:)), ...
        t ...
    );

    subplot(2,4,1); imshow(Y(:,:,1), []); title('X')
    subplot(2,4,2); imshow(Uc*(T(:,:,1) + M)*Ur', []); title('Rec. L');
    subplot(2,4,3); imshow(E(:,:,1), []); title('Rec. E');
    subplot(2,4,4); meshplot(Uc,1,1); title('Uc');
    subplot(2,4,5); meshplot(Ur,1,1); title('Ur');
    subplot(2,4,6); imshow(O(:,:,1), []); title('O')
    subplot(2,4,7); stemplot(T(:,:,1),1,1); title('T_1');
    subplot(2,4,8); meshplot(T(:,:,1),1,1); title('T_1');
    drawnow;
end

Data.T = T;
Data.Uc = Uc;
Data.Ur = Ur;
Data.E = E;
Info.gamma = gamma;
Info.alpha = alpha;

end

%% Tools that may be useful

function [ M ] = get_traces(A)
    r = size(A, 3);
    M = zeros(r);
    
    for i=1:r
        M(i,i) = trace(A(:,:,i));
    end
end

function ppup(Q, t)
    fprintf('Updated %s in %fs\n', Q, t);
end

% function print_condition_number(M)
%     [~, S, ~] = svd(M);
%     S = diag(S);
%     fprintf('Condition number of %s : %f\n', inputname(1), max(S)/min(S));
% end

function print_condition_number(M)
%     [~, S, ~] = svd(M);
%     S = diag(S);
    fprintf('Condition number of %s : %f\n', inputname(1), cond(M));
end
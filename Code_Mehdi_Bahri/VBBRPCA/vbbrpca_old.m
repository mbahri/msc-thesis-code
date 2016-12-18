function [ Uc, Ur, T, E ] = vbbrpca( Y, T_in, E_in, Uc_in, Ur_in, r, MAXITER)
%VBBRPCA VBBRPCA with Uc and Ur column independent

addpath Bayes_tools/

% MAXITER = 1;

% Constants about the data
[n, m] = size(Y(:,:,1));
Nobs = size(Y, 3);
% M = mean(T_in, 3);
M = zeros(r);

% Initialization of gamma, alpha, beta
a = 1e-6;
b = 1e-6;

Y2sum = sum(Y(:).^2);
scale2 = Y2sum / (Nobs*m*n);
scale = sqrt(scale2);

% gamma = gamrnd(a, b, 1, r);
gamma = scale*ones(1,r);

% gamma = scale*ones(r,1);        
beta = 1./scale2;
alpha = ones(n,m)*scale;

% Work variables
T = T_in;
E = E_in;
Uc = Uc_in;
Ur = Ur_in;

% Temp vars
mu_T = zeros(r^2, Nobs);

Sigma_Uc_i = scale*ones(n, n, r);
Sigma_Ur_i = scale*ones(m, m, r);

% Main loop
for it = 1:MAXITER,    
    % We use the expectations of Uc'Uc and Ur'Ur a lot so we precompute
    % them
%     UcUc = Uc'*Uc + sum(Sigma_Uc_i, 3);
%     UrUr = Ur'*Ur + sum(Sigma_Ur_i, 3);
    UcUc = Uc'*Uc + get_traces(Sigma_Uc_i);
    UrUr = Ur'*Ur + get_traces(Sigma_Ur_i);
    
    % Update of En and Tn
    beta_alpha_ratio = ( beta ./ (beta + alpha) );
    Sigma_T = inv( (eye(r^2) + beta*kron(UrUr, UcUc)) );

    % Temporary variables - preallocation
    for k=1:Nobs
        % Every E(i,j,k) can be computed at once by this vectorized
        % expression
        E(:,:,k) = beta_alpha_ratio .* ( Y(:,:,k) - Uc*T(:,:,k)*Ur' );
        
       % Compute mu1'
       mu_T(:,k) = reshape( M + beta*Uc'*(Y(:,:,k) - E(:,:,k))*Ur , r^2, 1); 
       mu_n = Sigma_T *  mu_T(:, k);
       T(:,:,k) = reshape( mu_n, r, r );
    end
    
    % Update Uc
    for i = 1:r
        Si = extract_row(i, r);
        Sigma_ni = Si*Sigma_T*Si';
        
        Sum_wni_wni = 0;
        mu_temp = zeros(n, 1);
        
        % Compute the sum over n
        for k=1:Nobs
            Sum_k_not_i = zeros(n);
            
            mu_ni = Si*mu_T(:, k);
            w_ni = Ur*mu_ni;
            Sum_wni_wni = Sum_wni_wni + trace( UrUr * (Sigma_ni + mu_ni*mu_ni') );
            
            for j=1:r
                if j ~= i
                    Sum_k_not_i = Sum_k_not_i + Uc(:,j) * (Ur * extract_row(j, r) * mu_T(:, k))' ;
                end
            end
            
            mu_temp = mu_temp + ( Y(:,:,k) - E(:,:,k) - Sum_k_not_i ) * w_ni;
        end
        
        Sigma_temp = (gamma(i) + beta * Sum_wni_wni)^-1 * eye(n);
        Sigma_Uc_i(:,:,i) = Sigma_temp;
        Uc(:,i) = Sigma_temp * (beta * mu_temp);
    end
    
    % Update Ur
    for i = 1:r
        Si = extract_column(i, r);
        Sigma_ni = Si*Sigma_T*Si';
        
        Sum_wni_wni = 0;
        mu_temp = zeros(m, 1);
        
        % Compute the sum over n
        for k=1:Nobs
            Sum_k_not_i = zeros(m);
            
            mu_ni = Si*mu_T(:, k);
            w_ni = Uc*mu_ni;
            Sum_wni_wni = Sum_wni_wni + trace( UcUc * (Sigma_ni + mu_ni*mu_ni') );
            
            for j=1:r
                if j ~= i
                    Sum_k_not_i = Sum_k_not_i + Ur(:,j) * (Uc * extract_column(j, r) * mu_T(:, k))' ;
                end
            end
            
            mu_temp = mu_temp + ( Y(:,:,k) - E(:,:,k) - Sum_k_not_i )' * w_ni;
        end
        
        Sigma_temp = (gamma(i) + beta * Sum_wni_wni)^-1 * eye(n);
        Sigma_Ur_i(:,:,i) = Sigma_temp;
        Ur(:,i) = Sigma_temp * (beta * mu_temp);
    end
    
    % Update alpha
    alpha = Nobs ./ (Nobs * beta_alpha_ratio + sum(E, 3).^2);
    
    % Update beta    
    sum_of_norms = 0;
    for k=1:Nobs
        sum_of_norms = sum_of_norms + ...
            expectation_norm(Y(:,:,k), Uc, T(:,:,k), Ur, E(:,:,k), UcUc, UrUr, Sigma_T, V, beta_alpha_ratio);
    end
%     expectation_norm( Y, Uc, T, Ur, E, UcUc, UrUr, Sigma_T, Var_E )
    
    beta = Nobs*m*n / sum_of_norms;
    
    % Update gamma
    gamma_b = 2*b ./ (2 + b * diag(UcUc + UrUr)); 
    gamma = (a + (m + n)/2) ./ gamma_b;
    
    stemplot(Uc,1,1);
    waitforbuttonpress;
    stemplot(Ur,1,1);
    waitforbuttonpress;
    fprintf('Continuing\n');
end

end

function [ M ] = get_traces(A)
    r = size(A, 3);
    M = zeros(r);
    
    for i=1:r
        M(i,i) = trace(A(:,:,i));
    end
end

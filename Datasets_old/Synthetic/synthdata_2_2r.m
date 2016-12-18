% clear;

r_c = 23;
r_r = 31;
n = 1000;
m = 1000;

Nobs = 100;
X = zeros(n, m, Nobs);
O = zeros(n, m, Nobs);
E_ori = zeros(n, m, Nobs);

T_ori = zeros(r_c, r_r, Nobs);
Uc_ori = rand(n, r_c); 
Ur_ori = rand(m, r_r);

randn('state',1212412414424234324);
rand('state',1212412414424234324);

% Uc_orig = randn(n,r);
% Ur_orig = randn(r,m);

for i=1:Nobs
    T = randn(r_c, r_r);
    T_ori(:,:,i) = T;
    D0 = Uc_ori*T*Ur_ori' / sqrt(r_c * r_r);
    
    % Save original D
    O(:,:,i) = D0;

    E0 = sign(randn(n,m));
    inds = rand(n,m)<0.7;
    E0(inds) = 0;
    
    % Save original E
    E_ori(:,:,i) = E0;

    X(:,:,i) = D0 + E0;
end
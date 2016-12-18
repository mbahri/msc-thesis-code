function [ Ori, X ] = synthetic_data_30( m, n, Nobs, r_c, r_r )
%SYNTHETIC_DATA_30 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5
    r_r = r_c;
end

X = zeros(n, m, Nobs);
Ori.O = zeros(n, m, Nobs);
Ori.E = zeros(n, m, Nobs);

Ori.T = zeros(r_c, r_r, Nobs);

randn('state',1212412414424234324);
rand('state',1212412414424234324);

Ori.Uc = rand(n, r_c); 
Ori.Ur = rand(m, r_r);

for i=1:Nobs
    T = randn(r_c, r_r);
    Ori.T(:,:,i) = T;
    D0 = Ori.Uc*T*Ori.Ur' / sqrt(r_c * r_r);
    
    % Save original D
    Ori.O(:,:,i) = D0;

    E0 = sign(randn(n,m));
    inds = rand(n,m)<0.7;
    E0(inds) = 0;
    
    % Save original E
    Ori.E(:,:,i) = E0;

    X(:,:,i) = D0 + E0;
end


end


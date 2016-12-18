clear;

r = 20;
n = 1000;
m = 1000;

Nobs = 10;
X = zeros(n, m, Nobs);
O = zeros(n, m, Nobs);
% E_ori = zeros(n, m, Nobs);

T_ori = rand(r, r, Nobs);
Uc_ori = rand(n, r); 
Ur_ori = rand(m, r);

for i=1:Nobs
    A_ymp = Uc_ori * T_ori(:,:,i) * Ur_ori';
    O(:,:,i) = A_ymp;
    
    Mask = randi([0,100], n, m);
    whites = Mask >= 0 & Mask <= 5;
    blacks = Mask <= 100 & Mask >= 95;
    A_ymp(whites) = 0;
    A_ymp(blacks) = 255;
    
%     E_tmp = ones(n, m);
%     E_tmp(whites) = 0;
%     E_tmp(blacks) = 255;
%     
%     E_ori(:,:,i) = E_tmp;
    X(:,:,i) = A_ymp;
end
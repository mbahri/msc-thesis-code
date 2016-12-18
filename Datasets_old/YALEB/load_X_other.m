load('yaleb10.mat');
X = X(:,:,1:64);

% X = X - repmat(mean(X,3), 1, 1, 64);

Nobs = size(X, 3);
F = 4;

[h, w] = size(X(:,:,1));
h1 = F*h; w1 = F*w;

XX = X;
% X = zeros(h1, w1, Nobs);
O = zeros(h1, w1, Nobs);

for i=1:Nobs
    A = imresize(XX(:,:,i), F);
    O(:,:,i) = A;
    
%     Mask = randi([0,100], h1, w1);
%     
%     whites = Mask >= 0 & Mask <= 29;
%     blacks = Mask <= 100 & Mask >= 71;
%     
%     A(whites) = -255;
%     A(blacks) = 255;  
%     
%     X(:,:,i) = A;
end
X = imnoise(O,'salt & pepper', 0.1);
Z = double(tenmat(tensor(X), 3))';

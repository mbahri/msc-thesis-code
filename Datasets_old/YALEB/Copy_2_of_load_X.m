load('YaleB_32x32.mat');
% load('5Train/1.mat');

% trainIdx = 1:64;

% feat = fea(trainIdx, :);
feat = fea;
Nobs = size(feat,1);
X = cell(1, Nobs);
O = cell(1, Nobs);
whites = cell(1, Nobs);
blacks = cell(1, Nobs);

for i=1:Nobs
    A = reshape(feat(i,:), 32, 32);
    A = imresize(A, 4);
    A = A ./ max(abs(A(:)));

    O{i} = A;

    Mask = randi([0,100], 128, 128);
%     Mask = randi([0,100], 32, 32);
    whites{i} = Mask >= 0 & Mask <= 29;
    blacks{i} = Mask <= 100 & Mask >= 71;
    A(whites{i}) = 0;
    A(blacks{i}) = 255;
    
    X{i} = A;
end

X=cat(3, X{:});
O=cat(3, O{:});
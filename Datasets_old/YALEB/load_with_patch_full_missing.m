load('YaleB_32x32.mat');
load('5Train/1.mat');


% feat = fea(1:45, :);
feat = fea;

N = size(feat, 1);
X = cell(1, N);
O = cell(1, N);
whites = cell(1, N);
blacks = cell(1, N);

for i=1:N
    A = reshape(feat(i,:), 32, 32);
    A = imresize(A, 4);

    O{i} = A;

%     Mask = randi([0,100], 50, 50);
%     Mask = randi([0,100], 32, 32);
%     A(40:89,40:89) = 0;
%     A(40:89,40:89) = max(A(40:89,40:89) - 255*(Mask >= 80 & Mask <= 100), 0);
   
    X{i} = A;
end

A = X{1};
A(60:69,60:69) = 0;
X{1} = A;

X=cat(3, X{:});
O=cat(3, O{:});
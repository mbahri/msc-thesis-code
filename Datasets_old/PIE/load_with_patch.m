load('fea64');
load('Train5_64');

% trainidx = Train5_64(1, :);
% feat = fea64(trainidx, :);
feat = fea64;

N = size(feat, 1);
X = cell(1, N);
O = cell(1, N);
whites = cell(1, N);
blacks = cell(1, N);

for i=1:N
    A = reshape(feat(i,:), 64, 64);
%     A = imresize(A, 2);

    O{i} = A;

    Mask = randi([0,100], 20, 20);
%     Mask = randi([0,100], 32, 32);
    A(23:42,23:42) = min(A(23:42,23:42) + 255*(Mask >= 0 & Mask <= 10), 255);
%     A(23:42,23:42) = max(A(23:42,23:42) - 255*(Mask >= 90 & Mask <= 100), 0);
   
    X{i} = A;
end

X=cat(3, X{:});
O=cat(3, O{:});
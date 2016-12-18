load('fea64');
load('Train5_64');

% trainidx = Train5_64(1, :);
trainidx = 1:49;
feat = fea64(trainidx, :);

N = size(trainidx, 2);

X = cell(1, N);
O = cell(1, N);
whites = cell(1, N);
blacks = cell(1, N);

for i=1:N
    A = reshape(feat(i,:), 64, 64);
    A = imresize(A, 2);

    O{i} = A;

    Mask = randi([0,100], 128, 128);
%     Mask = randi([0,100], 32, 32);
    whites{i} = Mask >= 0 & Mask <= 4;
    blacks{i} = Mask <= 100 & Mask >= 96;
    A(whites{i}) = 0;
    A(blacks{i}) = 255;
    
    X{i} = A;
end

X=cat(3, X{:});
O=cat(3, O{:});
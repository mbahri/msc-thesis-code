load('yaleb10.mat');

orig = cell(1, 640);
XX = cell(1, 640);
whites = cell(1, 640);
blacks = cell(1, 640);

for i=1:640
    A = reshape(X(:,i), imsize);
    A = imresize(A, 2);
    orig{i} = A;
    Mask = randi([0,100], 96, 84);
%     Mask = randi([0,100], 32, 32);
    whites{i} = Mask == 0 | Mask == 1 | Mask == 2 | Mask == 3 | Mask == 4;
    blacks{i} = Mask == 100 | Mask == 99 | Mask == 98 | Mask == 97 | Mask == 96;
    A(whites{i}) = 0;
    A(blacks{i}) = 255;
    XX{i} = A;
end
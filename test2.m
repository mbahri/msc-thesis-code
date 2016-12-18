clear

init
synthdata_2
% load_X
% E_ori = zeros(128,128,64);

taille = 1000;

YY = zeros(taille, taille*taille);
%OO = zeros(taille, taille*taille);
%EE = zeros(taille, taille*taille);

for i=1:100
    YY(i,:) = to_vec(X(:,:,i));
%    OO(i,:) = to_vec(O(:,:,i));
%    EE(i,:) = to_vec(E_ori(:,:,i));
end

fprintf('\n****** STARTING ******\n');

%options.X_true = OO';
%options.E_true = EE';
options.verbose = 1;

[Xrec, Arec, Brec, Erec] = VBRPCA(YY', options);

fprintf('\n****** DONE TRANSPOSE ******\n');

%options2.X_true = OO;
%options2.E_true = EE;
options2.verbose = 1;

[Xrec2, Arec2, Brec2, Erec2] = VBRPCA(YY, options2);

fprintf('\n****** DONE NO TRANSPOSE ******\n');

options1.X_true = O(:,:,1);
options1.E_true = E_ori(:,:,1);
options1.verbose = 1;

[X1, A1, B1, E1] = VBRPCA(X(:,:,1), options1);

fprintf('\n****** DONE FIRST MATRIX ******\n');

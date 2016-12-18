function [ output_args ] = show4a( n, X, Ori, Uc, Ur, T, E )
%SHOW4a n, X, Ori, Uc, Ur, T, E

subplot(1,4,1), imshow(X(:,:,n), [])
subplot(1,4,2), imshow(Ori(:,:,n), [])
subplot(1,4,3), imshow(Uc*T(:,:,n)*Ur', [])
subplot(1,4,4), imshow(E(:,:,n), [])
end


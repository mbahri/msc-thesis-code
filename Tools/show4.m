function [ output_args ] = show4( n, X, Ori, Uc, Ur, T, E )
%SHOW4 n, X, Ori, Uc, Ur, T, E
K = Ori(:,:,n);
a = min(K(:));
b = max(K(:));
subplot(1,4,1), imshow(X(:,:,n), [a,b])
subplot(1,4,2), imshow(Ori(:,:,n), [a,b])
subplot(1,4,3), imshow(Uc*T(:,:,n)*Ur', [a,b])
subplot(1,4,4), imshow(E(:,:,n), [a,b])
end


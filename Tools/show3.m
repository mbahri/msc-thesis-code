function [ output_args ] = show3( n, X, Uc, Ur, T, E )
subplot(1,3,1), ddisp(X(:,:,n))
subplot(1,3,2), ddisp(Uc*T(:,:,n)*Ur')
subplot(1,3,3), ddisp(E(:,:,n))
end


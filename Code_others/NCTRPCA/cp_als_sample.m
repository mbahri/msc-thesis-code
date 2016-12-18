function [w,A,B,C] = cp_als_sample(X1,X2,X3,r)
A = rand(size(X1,1),r);
B = rand(size(X2,1),r);
C = rand(size(X3,1),r);
w = ones(r,1);
w_prev = zeros(r,1);
CP_TOL = 1e-17;
MAX_ITER = 100;
t = 0
while norm(w-w_prev)>CP_TOL && t<MAX_ITER
    t = t+1
    w_prev = w;
    A = X1*((X3'*C).*(X2'*B))/size(X1,2);
    A = A./repmat(sqrt(sum(A.*A,1)),size(X1,1),1);
    B = X2*((X1'*A).*(X3'*C))/size(X1,2);
    B = B./repmat(sqrt(sum(B.*B,1)),size(X2,1),1);
    C = X3*((X1'*A).*(X2'*B))/size(X1,2);
    w = sqrt(sum(C.*C,1));
    w = w';
    C = C./repmat(sqrt(sum(C.*C,1)),size(X3,1),1);
end
end

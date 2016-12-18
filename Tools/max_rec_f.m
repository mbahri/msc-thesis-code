function [ m ] = max_rec_f( Uc, T, Ur, E, X )
%MAX_REC_F Max reconstruction error for the full observation
N = size(T,3);
m = -Inf;
for i=1:N
    rec = norm( Uc*T(:,:,i)*Ur' + E(:,:,i) - X(:,:,i), 'fro' ) / norm( X(:,:,i), 'fro' );
    if rec > m
        m = rec;
    end
end
end


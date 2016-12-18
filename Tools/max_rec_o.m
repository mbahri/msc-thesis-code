function [ m ] = max_rec_o( Uc, T, Ur, O )
%MAX_REC_O Max reconstruction error for the low rank term
N = size(T,3);
m = -Inf;
for i=1:N
    rec = norm( Uc*T(:,:,i)*Ur' - O(:,:,i), 'fro' ) / norm( O(:,:,i), 'fro' );
    if rec > m
        m = rec;
    end
end
end


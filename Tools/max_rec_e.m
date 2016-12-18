function [ m ] = max_rec_e( E, E_ori )
%MAX_REC_E Max reconstruction error for the sparse error
N = size(E,3);
m = -Inf;
for i=1:N
    rec = norm( E(:,:,i) - E_ori(:,:,i), 'fro' ) / norm( E_ori(:,:,i), 'fro' );
    if rec > m
        m = rec;
    end
end
end


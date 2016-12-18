function [ V ] = get_V(Sigma_T, r)

V = col2im(Sigma_T, [r, r], [r^2, r^2], 'distinct');
% V = blocktranspose(V, r^2, r^2, r, r);

end

function [ V ] = get_V_p(Sigma_T, r1, r2)

V = col2im(Sigma_T, [r1, r2], [r1*r2, r1*r2], 'distinct');
% V = blocktranspose(V, r^2, r^2, r, r);

end

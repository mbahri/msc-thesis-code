tT = tensor(T_ori);

res1 = zeros(n,m,Nobs);
for i=1:Nobs
    res1(:,:,i) = Uc_ori*T_ori(:,:,i)*Ur_ori';
end

res2 = ttm(tT, {Uc_ori, Ur_ori}, [1 2]);
% res3 = ttensor(T_ori, Uc_ori, Ur_ori, eye(100));
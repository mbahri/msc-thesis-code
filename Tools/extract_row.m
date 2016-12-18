function [ S ] = extract_row( i, r )
%EXTRACT_ROW Returns the matrix that extract row i from vec(Tn)

S = zeros(r, r^2);
for k=1:r
    S(k, i + (k-1)*r) = 1;
end

end


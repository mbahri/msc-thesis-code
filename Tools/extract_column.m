function [ S ] = extract_column( i, r )
%EXTRACT_COLUMN Returns the matrix that extract col i from vec(Tn)

S = [ repmat(zeros(r), 1, i-1) eye(r) repmat(zeros(r), 1, r-i) ];

end


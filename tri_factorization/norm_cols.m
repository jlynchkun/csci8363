function [A,m] = norm_cols(A,ell)

for i = 1:size(A,2)
    m(i) = norm(A(:,i),ell);
    if m(i) > 0
        A(:,i) = A(:,i)/m(i);
    end
end
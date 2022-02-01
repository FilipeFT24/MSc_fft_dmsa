function [X] = gaussjordanReduction(M)
[m,n] = size(M);
j     = 1;
k     = 1;
while (j <= m && k <= n)
    [r,c] = find(M(j:n,k:m),1);
    if ~isempty(r) && ~isempty(c)
        if c > 1
            k = k+c-1;
        end
        if r > 1
            M([r,j],:) = M([j,r],:);
        end
        M(j,:)    = M(j,:)/M(j,k);
        rows      = [1:j-1,j-1:m];
        M(rows,:) = M(rows,:)-M(rows,k).* M(j,:);
        
    end
    j = j+1;
    k = k+1;
end
X = M;
end
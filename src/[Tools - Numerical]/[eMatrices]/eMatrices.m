function [Inverse] = eMatrices(A)
%[E,P,S,Inverse] = eMatrices(A)
%
% This program find the elimination matrices of which will reduce any A
% matrix using the Gauss-Jordan method without partial pivoting. It has
% three outputs. It then goes on to find the inverse of the A matrix using
% the elimination matrices
%
% E = elimination matrices which reduce A to U
% P = elimination matrices which reduce U to a diagonal matrix
% S = scaling matrices which reduce the diagonal matrix to an identity
% matrix

% Written by Olawale B. Akinwale
% Department of Electronic and Electrical Engineering,
% Obafemi Awolowo University
% 2010-06-13

[m,n] = size(A);
if (m ~= n)
    return;
end
x = 1; % iterator for elimination matrix 3rd dimension
for j = 1:1:n-1
    for i = j+1:1:m
        E(:,:,x) = eye(n);
        E(i,j,x) = -A(i,j) / A(j,j);
        A = E(:,:,x) * A;
        x = x + 1;
    end
end
x = 1;
for j = n:-1:2
    for i = j-1:-1:1
        P(:,:,x) = eye(n);
        P(i,j,x) = -A(i,j) / A(j,j);
        A = P(:,:,x) * A;
        x = x + 1;
    end
end
for i = 1:1:n
    S(:,:,i) = eye(n);
    S(i,i,i) = 1 / A(i,i);
    A = S(:,:,i) * A;
end

ProdE = 1;
[~,~,c] = size(E);
for ii = c:-1:1
    ProdE = ProdE * E(:,:,ii);
end
ProdP = 1;
[~,~,c] = size(P);
for ii = c:-1:1
    ProdP = ProdP * P(:,:,ii);
end
ProdS = 1;
[~,~,d] = size(S);
for jj = d:-1:1
    ProdS = ProdS * S(:,:,jj);
end
Inverse = ProdS*ProdP*ProdE;
end
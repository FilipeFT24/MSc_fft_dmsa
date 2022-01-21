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

tic;

[m,n]    =  size(A);
x        =  1;          % iterator for elimination matrix 3rd dimension
j        =  1:1:n-1;
i        =  j+1:1:m;
E(:,:,x) =  eye(n);
E(i,j,x) = -A(i,j)/A(j,j);
A        =  E(:,:,x)*A;
x        =  x+1;

x        =  1;
j        =  n:-1:2;
i        =  j-1:-1:1;
P(:,:,x) =  eye(n);
P(i,j,x) = -A(i,j)/A(j,j);
A        =  P(:,:,x)*A;
x        =  x+1;

for i = 1:1:n
    S(:,:,i) = eye(n);
    S(i,i,i) = 1 / A(i,i);
    A = S(:,:,i) * A;
end

toc;




ProdE    = 1;
[~,~,c]  = size(E);
ii       = c:-1:1;
ProdE    = ProdE*E(:,:,ii);

ProdP    = 1;
[~,~,c]  = size(P);
ii       = c:-1:1;
ProdP    = ProdP*P(:,:,ii);

ProdS    = 1;
[~,~,d]  = size(S);
jj       = d:-1:1;
ProdS    = ProdS*S(:,:,jj);
Inverse  = ProdS*ProdP*ProdE;
end
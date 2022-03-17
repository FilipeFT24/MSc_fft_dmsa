function [Inverse] = eMatrices(A)
[m,n] = size(A);
 
I = eye(m);
 
B = [A,I];
 
[m,n] = size(B);
 
% All row operations
 
% Downwards
 
for i=1:m 
 
B(i,:) = B(i,:)/B(i,i);
 
for j=(i+1):m 
 
B(j,:) = B(j,:)+(-1)*B(j,i)*B(i,:);
 
end 
 
end 
 
%Upwards
 
for i=m:-1:1
 
for j=(i-1):-1:1
 
B(j,:) = B(j,:) + (-1)*B(j,i)*B(i,:);
 
end 
 
end 
 
 

%Inverse
 
Inverse = B(:,(n/2)+1:n);
end
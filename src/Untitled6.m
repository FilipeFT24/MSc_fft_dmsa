clear; clc;

A =[1 1 2 2 3 4 5];

m = find(accumarray(A.',ones(size(A))) > 1)
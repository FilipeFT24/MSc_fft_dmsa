clear; clc;

A = [13,24;24,25;14,25;13,14]
B = [13,14;13,24;24,25;14,25]


[~,i_AS] = sortrows(A); 
[~,j]    = sortrows(i_AS);
[~,i_BS] = sortrows(B);


ret = i_BS(j)

B(i_BS(j),:)



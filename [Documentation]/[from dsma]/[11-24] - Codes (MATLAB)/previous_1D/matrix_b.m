function [b]=matrix_b
% Fun��o que obtem a matriz preenche a matrix b com os valores obtidos pelo
% 'integral de gauss'



%% Declara��o das Variaveis Globais %%
%
global L Lref cell_num face_num;
global x faces cell_vol lap_phi;
%
%% Constru��o da Matriz %%

% A=zeros(cell_num);
b=sparse([]);

for i=1:cell_num
        b(i)=lap_phi(i)*cell_vol(i);
end
%
end
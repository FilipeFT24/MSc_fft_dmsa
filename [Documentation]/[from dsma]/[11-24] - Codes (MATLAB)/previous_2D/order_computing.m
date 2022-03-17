clc


%% ORDEM
tamanho=size(vector_save4_1);
for iter=2:tamanho(1)

tamanho_malha1= vector_save8_1(iter-1,1)
tamanho_malha2= vector_save8_1(iter,1)

norm1_1=vector_save8_1(iter-1,2);
norm1_2=vector_save8_1(iter,2);

normax_1=vector_save8_1(iter-1,3);
normax_2=vector_save8_1(iter,3);

ordem_convergencia8(iter)= log(abs((norm1_2/norm1_1)))/log(abs((tamanho_malha2/tamanho_malha1)));
ordem_convergencia8_max(iter)= log(abs((normax_2/normax_1)))/log(abs((tamanho_malha2/tamanho_malha1)));
end
%% RACIO
% 
% tamanho=size(vector_save4_1);
% 
% for i=1:tamanho(1)
%    racio_med(i)= vector_save4_2(i,2)/vector_save4_1(i,2);
%    racio_max(i)= vector_save4_2(i,3)/vector_save4_1(i,3);
% end







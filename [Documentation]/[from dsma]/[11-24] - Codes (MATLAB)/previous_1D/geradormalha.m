function [x,faces,cell_vol]=geradormalha(type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   07 de Junho de 2016                                             %
%                                                                                                   %
% Fun��o que cria uma malha uniforme e n�o uniforme                                                 %
% Malha n�o uniforme � criada com base na malha uniforme com deslocamento das faces para a esquerda %
% ou para a direita de forma aleat�ria bem como a distancia que a face � deslocada, sempre tendo por%
% base o comprimento de refencia da c�lula.                                                         %
%                                                                                                   %
% x         - Coordenadas dos centroides das c�lulas;                                               %
% faces     - Coordenadas dos centroides das faces;                                                 %
% cell_vol  - Volume das c�lulas.                                                                   %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Declara��o das Variaveis Globais %%
%
global L Lref cell_num face_num
%
%% Contru��o da Malha Uniforme %%
%
if type==true
    % Coordenadas das Faces das C�lulas %
    for i=1:face_num
        faces(i)=(i-1)*Lref;
    end
    %
    % Coordenadas dos Centroides das C�lulas %
    for i=1:cell_num
        if i==1
            x(i)=0+(faces(i+1)-0)/2;
        else
            x(i)=faces(i)+(faces(i+1)-faces(i))/2;
        end
    end
    %
%% Constu��o da Malha N�o Uniforme %%
%
else
    % Variaveis Auxiliares para a constru��o da malha n�o uniforme %
    % aux1 -> Coordenadas das Faces de uma malha regular %
    for i=1:face_num
        aux1(i)=(i-1)*Lref;
    end
    %
    % aux2 -> Vector Random de inteiros entre 1 e 2 para determinar o lado para que a face se 
    % desloca, 1-> Face move para a esquerda, 2-> Face move para a direita.
    aux2=randi(2,face_num,1);
    %
    % aux3 -> Vector Random de inteiros entre 1 e 2 para determinar o tamanho do deslocamento da
    % face, 1-> Face desloca-se 0.10*Lref, 2 -> Face desloca-se 0.15*Lref
    aux3=randi(2,face_num,1);
    %
    % Determina��o das coordenadas das faces %
    for i=1:face_num
        if i==1
            faces(i)=0;
        elseif i==face_num
            faces(i)=L;
        else
            if aux2(i)==1
                if aux3(i)==1
                    faces(i)=aux1(i)-0.1*Lref;
                else
                    faces(i)=aux1(i)-0.15*Lref;
                end
            else
                if aux3(i)==1
                    faces(i)=aux1(i)+0.1*Lref;
                else
                    faces(i)=aux1(i)+0.15*Lref;
                end
            end
        end
    end
    %
    % Determina��o dos Centroides das C�lulas %
    for i=1:cell_num
        x(i)=(faces(i+1)-faces(i))/2+faces(i);
    end  
end
%
%% Calcula o Volume das C�lulas %%
%
for i=1:cell_num
   cell_vol(i)=faces(i+1)-faces(i); 
end
%
%% Escreve para um ficheiro as coordenadas das celulas e das faces %%
%
% filetitle=['cells',num2str(cell_num),'.dat'];
% fid=fopen(filetitle,'w');
% for i=1:cell_num
%     fprintf(fid,'%d %E\n',i,x(i));
% end
% fclose(fid);
% %
% filetitle=['faces',num2str(cell_num),'.dat'];
% fid=fopen(filetitle,'w');
% for i=1:cell_num+1
%     fprintf(fid,'%d %E\n',i,faces(i));
% end
% fclose(fid);
% %
end
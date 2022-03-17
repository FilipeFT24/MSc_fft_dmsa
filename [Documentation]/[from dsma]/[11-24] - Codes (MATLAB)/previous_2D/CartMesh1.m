function [verts,cells,faces,cell_vert,cell_face,cell_vol,face_area,face_vert,cell_norm,cell_num,vert_num,face_num,Lref]=CartMesh1(uniforme)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos
%                                  16 de Outubro de 2016
%
% Função que cria uma malha Cartesiana
%
% verts     - Vetor que guarda as coordenadas dos vertices da malha         - verts(i)=(coordenada x, coordenada y);
% face_vert - Vetor que guarda os vertices das faces                        - face_vert(i)=(vertice 1, vertice 2);
% face_area - Vetor que guarda as areas dos vertices das faces              - face_area(i)=(area, area em x, area em y);
% faces     - Vetor que guarda as coordenadas do centroide da face          - faces(i)=(coordenada x, coordenada y)
% cell_vert - Vetor que guarada os vertices que compõem a celula            - cell_vert(i)=(vertice 1,vertice2,vertice4,vertice4)
% cell_face - Vetor que guarda as faces que compõem a celula                - cell_face(i)=(face 1,face 2,face 3,face 4)
% cell_vol  - Vetor que guarda o volume das células                         - cell_vol(i)=volume célula
% cells     - Vetor que guarda as coordendas do centroide da célula         - cells(i)=(coordenda x, coordenada y)
% cell_norm - Vetor que guarad as normais exteriores de cada face da célula - cell_norm(i)=(ne_x1,ne_y1,ne_x2,ne_y2,ne_x3,ne_y3,ne_x4,ne_y4,) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Declaração das Variaveis Globais %%
%
global L cell_side vert_side;
%
%% Determinação de Proprieades Básicas %%
%
cell_num=cell_side*cell_side;       % Numero Total de Células %
vert_num=vert_side*vert_side;       % Numero Total de Vertices %
face_num=cell_side*vert_side*2;     % Numero Total de Faces %
%
Lref=sqrt(L*L/cell_num);            % Comprimento de Referencia da Malha %
%
%% Contrução da Malha Uniforme %%
%
% Coordenadas Gerais dos Vertices %
if uniforme
    for i=1:vert_side
        aux(i)=(i-1)*Lref;
    end
else
    error('\n\nERRO: Malha Não Implementada\n\n');
end
%
%% Determinação dos Vertices %%
%
k=0;
%
% Numeração dos vertices segundo a vertical %
for i=1:vert_side
    for j=1:vert_side
        k=k+1;
        verts(k,1)=aux(i);          % Coordenada X do Vertice %
        verts(k,2)=aux(j);          % Coordenada Y do Vertice %
    end
end
%
%% Determinação das Faces %%
%
k=0;
%
% Faces Verticais %
for i=1:vert_side
    for j=1:vert_side-1
        k=k+1;
        %
        % Vertices da Face %
        face_vert(k,1)=vert_side*(i-1)+j;
        face_vert(k,2)=vert_side*(i-1)+j+1;
        %
        % Normal Exterior %
        ne(k,1)=1;
        ne(k,2)=0;
    end
end
%
% Faces Horizontais %
for i=1:vert_side-1
    for j=1:vert_side
        k=k+1;
        %
        % Vertices da Face %
        face_vert(k,1)=vert_side*(i-1)+j;
        face_vert(k,2)=vert_side*(i-1)+j+vert_side;
        %
        % Normal Exterior %
        ne(k,1)=0;
        ne(k,2)=1;
    end
end
%
% Area das Faces e Coordenadas do Centroides das Faces %
for i=1:face_num
    x1=verts(face_vert(i,1),1);
    x2=verts(face_vert(i,2),1);
    y1=verts(face_vert(i,1),2);
    y2=verts(face_vert(i,2),2);
    %
    % Area da Face %
    face_area(i,1)=sqrt((x2-x1)^2+(y2-y1)^2);
    face_area(i,2)=face_area(i,1)*ne(i,2);
    face_area(i,3)=face_area(i,1)*ne(i,1);
    %
    % Coordenadas do Centroide da Face %
    faces(i,1)=(x1+x2)/2;
    faces(i,2)=(y1+y2)/2;
end
%
%% Determinação das Celulas %%
%
k=0;
%
% Determinação dos Vertices da Celula %
for i=1:vert_side-1
    for j=1:vert_side-1
        k=k+1;
        cell_vert(k,1)=vert_side*(i-1)+j;
        cell_vert(k,2)=vert_side*(i-1)+j+1;
    end
end
k=0;
for i=2:vert_side
    for j=1:vert_side-1
        k=k+1;
        cell_vert(k,4)=vert_side*(i-1)+j;
        cell_vert(k,3)=vert_side*(i-1)+j+1;
    end
end
%
% Determinação das Faces das Células %
for i=1:cell_num
        k1=0;
        k2=0;
        k3=0;
        k4=0;
    for j=1:face_num
        if face_vert(j,1)==cell_vert(i,1) && face_vert(j,2)==cell_vert(i,2)
            cell_face(i,1)=j;
            k1=1;
        end
        if face_vert(j,1)==cell_vert(i,2) && face_vert(j,2)==cell_vert(i,3)
            cell_face(i,2)=j;
            k2=1;
        end
        if face_vert(j,2)==cell_vert(i,3) && face_vert(j,1)==cell_vert(i,4)
            cell_face(i,3)=j;
            k3=1;
        end
        if face_vert(j,2)==cell_vert(i,4) && face_vert(j,1)==cell_vert(i,1)
            cell_face(i,4)=j;
            k4=1;
        end  
    end
    if k1+k2+k3+k4~=4
        fprintf('\nErro Célula %d\t',i);
    end
end
%
% Determinação dos Centroides das Células %
for i=1:cell_num
    x1=verts(cell_vert(i,1),1);
    y1=verts(cell_vert(i,1),2);
    x2=verts(cell_vert(i,2),1);
    y2=verts(cell_vert(i,2),2);
    x3=verts(cell_vert(i,3),1);
    y3=verts(cell_vert(i,3),2);
    x4=verts(cell_vert(i,4),1);
    y4=verts(cell_vert(i,4),2);
    %
    % Coordenadas do Centroide da Celula %
    cells(i,1)=(x1+x2+x3+x4)/4;
    cells(i,2)=(y1+y2+y3+y4)/4;
    %
    % Volume da Celula %
    aux1=(x1*y2+x2*y3+x3*y4+x4*y1);
    aux2=(y1*x2+y2*x3+y3*x4+y4*x1);
    cell_vol(i)=abs(aux1-aux2)/2;
    %
    % Normais Exteriores das Faces da Célula %
    cell_norm(i,1)=(y1-y2)/sqrt((x2-x1)^2+(y2-y1)^2);
    cell_norm(i,2)=(x2-x1)/sqrt((x2-x1)^2+(y2-y1)^2);
    %
    cell_norm(i,3)=(y2-y3)/sqrt((x3-x2)^2+(y3-y2)^2);
    cell_norm(i,4)=(x3-x2)/sqrt((x3-x2)^2+(y3-y2)^2);
    %
    cell_norm(i,5)=(y3-y4)/sqrt((x3-x4)^2+(y3-y4)^2);
    cell_norm(i,6)=(x4-x3)/sqrt((x3-x4)^2+(y3-y4)^2);
    %
    cell_norm(i,7)=(y4-y1)/sqrt((x4-x1)^2+(y4-y1)^2);
    cell_norm(i,8)=(x1-x4)/sqrt((x4-x1)^2+(y4-y1)^2);
    %
end
%
end
function [face_cells,vert_cells,face_bound,cell_bound,cell_vert_num, cell_face_num,vert_cell_num,vert_face_num]=CartMesh2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                 16 de Outubro de 2016                                             %
%                                                                                                   %
% Fun��o que determina as propriedades da malha                                                     %
%                                                                                                   %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global L Lref cell_side vert_side cell_num face_num vert_num;
global verts cells faces;
global cell_verts cell_faces cell_vol cell_norm;
global face_vert face_area;
%
%% Determina��o das Faces de Fronteira do Dominio %%


%Alterar o tipo de fronteiras
for i=1:face_num
    if faces(i,1)==0 || faces(i,1)==L
        face_bound(i,1)=1;
        face_bound(i,2)=1;
%         face_bound(i,4)=1;
        if faces(i,1)==L
        face_bound(i,4)=0;
        else
        face_bound(i,4)=0;
        end
    end
    if faces(i,2)==0 || faces(i,2)==L
        face_bound(i,1)=1;
        face_bound(i,3)=1;
        face_bound(i,4)=0;

    end
end
%
%% Determina��o das C�lulas de Fronteira do Dominio %%
%
n=size(cell_faces,2);
m=size(cell_verts,2);
%
for i=1:cell_num
    for j=1:4
        aux=cell_faces(i,j);
        if face_bound(aux,2)==1
            cell_bound(i,1)=1;
            cell_bound(i,2)=1;
        end
        if face_bound(aux,3)==1
            cell_bound(i,1)=1;
            cell_bound(i,3)=1;
        end
    end
    %
    cell_face_num(i)=0;
    for j=1:n
        if cell_faces(i,j)~=0
            cell_face_num(i)=cell_face_num(i)+1;
        end
    end
    %
    cell_vert_num(i)=0;
    for j=1:m
        if cell_verts(i,j)~=0
            cell_vert_num(i)=cell_vert_num(i)+1;
        end
    end
    %
end
%
%% Determina��o das C�lulas que cont�m a Face %%
%
for i=1:face_num
    z=0;
    for j=1:cell_num
        for k=1:4
            if cell_faces(j,k)==i
                z=z+1;
                face_cells(i,z)=j;
            end
        end
    end
end
%
%% Determina��o das C�lulas que cont�m o Vertice %%
%
for i=1:vert_num
    z=0;
    vert_cell_num(i)=0;
    for j=1:cell_num
        for k=1:4
            if cell_verts(j,k)==i
                z=z+1;
                vert_cells(i,z)=j;
                vert_cell_num(i)=vert_cell_num(i)+1;
            end
        end
    end
    %
    vert_face_num(i)=0;
    for j=1:face_num
        for k=1:2
            if face_vert(j,k)==i
                vert_face_num(i)=vert_face_num(i)+1;
            end
        end
    end
end
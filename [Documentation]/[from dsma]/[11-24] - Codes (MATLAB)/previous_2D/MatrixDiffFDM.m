function [A,source,source_faces]=MatrixDiffFDM(dirichlet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   04 de Julho de 2016                                             %
%                                                                                                   %
% Função que implementa a Solução Numérica a partir do Método de Diferenças Finitas de 2ª Ordem     %
% O Método Explicito ainda não foi implementado, sendo que o método implicito recorre à inversão da % 
% matriz                                                                                            %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global TempoGlobal;
global fid;
global L Lref cell_side vert_side cell_num face_num vert_num;
global verts cells faces;
global cell_verts cell_faces cell_vol cell_norm cell_bound cell_vert_num cell_face_num;
global face_vert face_cells face_area face_bound;
global vert_cells vert_cell_num vert_face_num;
global phi lap_phi  phi_faces flux_phi_faces;
%
fprintf('\n\n\t\tConstrução da Matriz A\t\t\t Inicio ...');
tempo_A=cputime;
%
cell=0;
for i=1:cell_num
    if rem(i*10,cell_num)==0
        fprintf(' %d%% ',i*100/cell_num);
    end
    cell=i;
    %
    % Faces da Célula %
    face_west=cell_faces(i,1);
    face_north=cell_faces(i,2);
    face_east=cell_faces(i,3);
    face_south=cell_faces(i,4);
    %
    % Area das Faces %
    area_west=face_area(face_west,1);
    area_east=face_area(face_east,1);
    area_south=face_area(face_south,1);
    area_north=face_area(face_north,1);
    %
    % Distancias Entre Centroides %
    if face_bound(face_west,1)==1
        % Face de Fronteira Oeste x=0 %
        if dirichlet
            d_west=cells(cell,1)-faces(face_west,1);
            cell_west=0;
        else
            error('\n\n\tERRO: Método Não Implementado\n\n');
        end
    else
        if face_cells(face_west,2)== cell
            cell_west=face_cells(face_west,1);
        else
            error('\n\n\tERRO: Esta celula não pode estar à esquerda da celula em estudo \n\n');
            %cell_west=face_cells(face_west,2);
        end
        d_west=cells(cell,1)-cells(cell_west,1);
    end
    %
    if face_bound(face_east,1)==1
        % Face de Fronteira Este x=L %
        if dirichlet
            d_east=faces(face_east,1)-cells(cell,1);
            cell_east=0;
        else
            error('\n\n\tERRO: Método Não Implementado\n\n');
        end
    else
        if face_cells(face_east,1)==cell
            cell_east=face_cells(face_east,2);
        else
            error('\n\n\tERRO: Esta celula não pode estar à direita da celula em estudo \n\n');
            %cell_east=face_cells(face_east,2);
        end
        d_east=cells(cell_east,1)-cells(cell,1);
    end
    %
    if face_bound(face_north,1)==1
        % Face de Fronteira Norte y=L %
        if dirichlet
            d_north=faces(face_north,2)-cells(cell,2);
            cell_north=0;
        else
            error('\n\nERRO: Método Não Implementado\n\n');
        end
    else
        if face_cells(face_north,1)==cell
            cell_north=face_cells(face_north,2);
        else
            error('\n\n\tERRO: Esta célula não pode estar em cima da celula em estudo\n\n');
            %
        end
        d_north=cells(cell_north,2)-cells(cell,2);
    end
    %
    if face_bound(face_south,1)==1
        % Face de Fronteira Sul y=0 %
        if dirichlet
            d_south=cells(cell,2)-faces(face_south,2);
            cell_south=0;
        else
            error('\n\nERRO: Método Não Implementado\n\n');
        end
    else
        if face_cells(face_south,2)==cell
            cell_south=face_cells(face_south,1);
        else
            error('\n\nERRO: Esta célula não pode estar abaixo da celula em estudo\n\n');
            %
        end
        d_south=cells(cell,2)-cells(cell_south,2);
    end
    %
    % Determinação dos Termos Fonte tendo já em conta as celulas de fronteira %
    % Determinação da Matriz A %
    a_aux(cell,cell)=-(area_west/d_west+area_east/d_east+area_north/d_north+area_south/d_south);
    %
    if cell_west==0
        source_west=area_west/d_west*phi_faces(face_west);
    else
        source_west=0;
        a_aux(cell,cell-cell_side)=area_west/d_west;
    end
    %
    if cell_east==0
        source_west=area_east/d_east*phi_faces(face_east);
    else
        source_east=0;
        a_aux(cell,cell+cell_side)=area_east/d_east;
    end
    %
    if cell_north==0
        source_north=area_north/d_north*phi_faces(face_north);
    else
        source_north=0;
        a_aux(cell,cell+1)=area_north/d_north;
    end
    %
    if cell_south==0
        source_south=area_south/d_south*phi_faces(face_south);
    else
        source_south=0;
        a_aux(cell,cell-1)=area_south/d_south;
    end
    source_faces(cell)=source_west+source_east+source_north+source_south;
    source(cell)=lap_phi(cell)*cell_vol(cell)-source_faces(cell);
end
%
fprintf('... Fim\n');
%
A=sparse(a_aux);
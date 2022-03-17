function [stencil_cells,stencil_faces,stencil_size]=Stencil(order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                  12 de Outubro de 2016                                            %
%                                  13 de Outubro de 2016                                            %
%                                                                                                   %
% Stencil de 6ª Ordem, é necessário expandir o stencil até aos 3 vizinhos de vertice de cada face,  %
% nas fronteiras é feito uma extensão do stencil de forma a que apenas seja extendido na direcção   %
% que o stencil não está completo, neste caso em cada direcção são necessários 6 pontos             %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global TempoGlobal fid;
global L Lref cell_side vert_side cell_num face_num vert_num;
global verts cells faces;
global cell_verts cell_faces cell_vol cell_norm cell_bound cell_vert_num cell_face_num extended_stencil;
global face_vert face_cells face_area face_bound;
global vert_cells vert_cell_num vert_face_num;
global phi lap_phi  phi_faces flux_phi_faces;
%
% Determinação do Tipo de Vizinhos %
% neib=0 2ª Ordem %
% neib=1 4ª Ordem %
% neib=2 6ª Ordem %
% neib=3 8ª Ordem %
%
neib=order/2-1;
%
face=0;
%
for i=1:face_num
    %
    face=i;
    stencil_cells_aux=0;
    z=0;
    vertices=0;
    %
    % Vizinhos de Vertice 1º Vizinhos %
    for j=1:2
        vert=face_vert(face,j);
        vertices(j)=vert;
        %
        % Adiciona as Celulas do 1º Vertice da Face %
        if j==1
            for k=1:vert_cell_num(vert)
                if vert_cells(vert,k)~=0
                    z=z+1;
                    stencil_cells_aux(z)=vert_cells(vert,k);
                end
            end
            n=size(stencil_cells_aux,2);
            %
            % Adiciona as Células do 2º Vertice da Face %
        else
            z=0;
            for k=1:vert_cell_num(vert)
                if vert_cells(vert,k)~=0
                    %
                    % Teste para não ser adicionadas Celulas Repetidas ao Stencil %
                    teste=0;
                    for q=1:n+z
                        if vert_cells(vert,k)==stencil_cells_aux(q)
                            teste=1;
                        end
                    end
                    %
                    % Adiciona a Celula do 2º Vertice da Face ao Stencil %
                    if teste==0
                        z=z+1;
                        stencil_cells_aux(n+z)=vert_cells(vert,k);
                    end
                end
            end
        end
    end
    %
    % Ciclo para os Restantes Vizinhos de Vertices da Face %
    %
    for j=1:neib
        %
        % Recolha dos Vertices do Stencil para Determinar os Proximos Vizinhos da Face %
        z=0;
        n=size(stencil_cells_aux,2);
        m=size(vertices,2);
        for k=1:n
            %
            % Celula Selecionada %
            cell=stencil_cells_aux(k);
            if cell~=0
                for q=1:cell_vert_num(cell)
                    %
                    % Teste para não adicionar vertices repetidos %
                    teste=0;
                    for p=1:m+z
                        if cell_verts(cell,q)==vertices(p)
                            teste=1;
                        end
                    end
                    %
                    % Adiciona o Vertice ao Vetor %
                    if teste==0
                        z=z+1;
                        vertices(m+z)=cell_verts(cell,q);
                    end
                end
            end
        end
        m=size(vertices,2);
        z=0;
        %
        % Vizinhos de Vertice J+1 Vizinhos %
        for k=1:m
            %
            % Vertice Selecionado %
            vert=vertices(k);
            %
            % Verifica se a Celula já está contina no Stencil %
            for q=1:vert_cell_num(vert)
                teste=0;
                if vert_cells(vert,q)~=0
                    for p=1:n+z
                        if vert_cells(vert,q)==stencil_cells_aux(p)
                            teste=1;
                        end
                    end
                    if teste==0
                        z=z+1;
                        stencil_cells_aux(n+z)=vert_cells(vert,q);
                    end
                end
            end
        end
    end
    %
    % Organizar o Stencil de Celulas por Ordem Decrescente %
    stencil_cells_aux=sort(stencil_cells_aux,'descend');
    stencil_size(face,2)=0;
    z=0;
    n=size(stencil_cells_aux,2);
    %
    for j=1:n
        if stencil_cells_aux(j)~=0
            z=z+1;
            stencil_cells(face,z)=stencil_cells_aux(j);
            stencil_size(face,2)=stencil_size(face,2)+1;
        end
    end
    %
    %% Adicionar as Faces de Fronteira ao Stencil %%
    %
    z=0;
    n=stencil_size(face,2);
    stencil_faces_aux=0;
    %
    for j=1:n
        %
        % Celula Selecionada %
        cell=stencil_cells(face,j);
        %
        % Verifica se a Celula Existe %
        if cell~=0
            for k=1:cell_face_num(cell)
                %
                % Seleciona uma Face da Celula %
                face_aux=cell_faces(cell,k);
                %
                % Verifica se a Face é de Fronteira %
                if face_bound(face_aux,1)==1
                    if z==0
                        z=z+1;
                        stencil_faces_aux(z)=face_aux;
                    else
                        teste=0;
                        for q=1:z
                            if stencil_faces_aux(q)==face_aux
                                teste=1;
                            end
                        end
                        %
                        if teste==0
                            z=z+1;
                            stencil_faces_aux(z)=face_aux;
                        end
                    end
                end
            end
        end
    end
    %
    % Organizar o Stencil das Faces %
    stencil_faces_aux=sort(stencil_faces_aux,'descend');
    stencil_size(face,3)=0;
    z=0;
    n=size(stencil_faces_aux,2);
    %
    for j=1:n
        if stencil_faces_aux(j)~=0
            z=z+1;
            stencil_faces(face,z)=stencil_faces_aux(j);
            stencil_size(face,3)=stencil_size(face,3)+1;
        end
    end
    %
    % Tamanho do Stencil %
    stencil_size(face,1)=stencil_size(face,2)+stencil_size(face,3);
    %
    %% Verificação da Necessidade de Expansao do Stencil %%
    %
    x_max=0;
    x_min=L;
    y_max=0;
    y_min=L;
    %
    area_stencil=0;
    %
    % Determina o Tamanho do Stencil a fim de determinar para que direcções é necessário expandir o
    % Stencil %
    %
    % Celulas do Stencil %
    for j=1:stencil_size(face,2)
        %
        % Celula Selecionada %
        cell=stencil_cells(face,j);
        %
        % Comprimento em X %
        if x_max<cells(cell,1)
            x_max=cells(cell,1);
        end
        if x_min>cells(cell,1)
            x_min=cells(cell,1);
        end
        %
        % Comprimento em Y %
        if y_max<cells(cell,2)
            y_max=cells(cell,2);
        end
        if y_min>cells(cell,2)
            y_min=cells(cell,2);
        end
        %
        area_stencil=area_stencil+cell_vol(cell);
    end
    %
    % Faces do Stencil %
    for j=1:stencil_size(face,3)
        %
        % Face Selecionada %
        face_aux=stencil_faces(face,j);
        %
        % Comprimento em X %
        if x_max<faces(face_aux,1)
            x_max=faces(face_aux,1);
        end
        if x_min>faces(face_aux,1)
            x_min=faces(face_aux,1);
        end
        %
        % Comprimento em Y %
        if y_max<faces(face_aux,2)
            y_max=faces(face_aux,2);
        end
        if y_min>faces(face_aux,2)
            y_min=faces(face_aux,2);
        end
    end
    %
    % Determinação do Comprimento de Referencia do Stencil - Importante para os casos onde a malha não é
    % uniforme ou a malha é não estruturada %
    L_stencil=sqrt(area_stencil/stencil_size(face,2));
    %
    % Determina em que direcção o Stencil deve crescer %
    L_x=abs(x_max-x_min)/L_stencil;
    L_y=abs(y_max-y_min)/L_stencil;
    %
    % Determina quantos pontos tem o Stencil em cada direcção %
    n_x=int8(L_x+1.1);
    n_y=int8(L_y+1.1);
    %
    incremento(face)=0;
    %
    %%
    % Normal Exterior da Face %
    
    
    vert=face_vert(face,:);
    nx=verts(vert(1),1)-verts(vert(2),1);
    ny=verts(vert(1),2)-verts(vert(2),2);
    
    if nx==0 && extended_stencil
        ordem_y=order+1;
        ordem_x=order;
    elseif ny==0 && extended_stencil
        ordem_x=order+1;
        ordem_y=order;
    else
        ordem_x=order;
        ordem_y=order;
    end
    
    
    %%
    
 
    
    while n_x<ordem_x | n_y<ordem_y
        %
        % Stencil tem de ser Extendido %
        incremento(face)=incremento(face)+1;
        vertices=0;
        z=0;
        %
        % Vector Auxiiar com todos os Vertices do Stencil %
        for j=1:stencil_size(face,2)
            %
            % Celula Selecionada %
            cell=stencil_cells(face,j);
            for k=1:cell_vert_num(cell)
                %
                % Vertice Selecionado %
                vert=cell_verts(cell,k);
                teste=0;
                for q=1:z
                    if vert==vertices(q)
                        teste=1;
                    end
                end
                %
                if teste==0
                    z=z+1;
                    vertices(z)=vert;
                end
            end
        end
        %
        stencil_cells_aux=0;
        stencil_cells_aux=stencil_cells(face,:);
        n=size(vertices,2);
        z=stencil_size(face,2);
        %
        % Adiciona Novas Celulas %
        for j=1:n
            %
            % Vertice Selecionado %
            vert=vertices(j);
            for k=1:vert_cell_num(vert)
                %
                % Celula Selecionada %
                cell=vert_cells(vert,k);
                %
                % Caso a Celula Exista %
                if cell~=0
                    %
                    % Verificação se a Celula é para adicionar ao Stencil %
                    teste=0;
                    %
                    for q=1:z
                        if cell==stencil_cells_aux(q)
                            teste=1;
                        end
                    end
                    %
                    if teste==1
                        continue
                    end
                    %
                    % Crescimento em XX %
                    if n_x<ordem_x
                        if n_y<order
                            teste=0;
                        else
                            if cells(cell,2)<y_min | cells(cell,2)>y_max
                                teste=1;
                            end
                        end
                    end
                    %
                    % Crescimento em YY %
                    if n_y<ordem_y
                        if n_x<order
                            teste=0;
                        else
                            if cells(cell,1)<x_min | cells(cell,1)>x_max
                                teste=1;
                            end
                        end
                    end
                    %
                    % Adiciona a Celula ao Stencil %
                    if teste==0
                        z=z+1;
                        stencil_cells_aux(z)=cell;
                    end
                end
            end
        end
        %
        % Organizar o Stencil %
        stencil_cells_aux=sort(stencil_cells_aux,'descend');
        stencil_size(face,2)=0;
        n=size(stencil_cells_aux,2);
        z=0;
        for j=1:n
            if stencil_cells_aux(j)~=0
                z=z+1;
                stencil_cells(face,z)=stencil_cells_aux(j);
                stencil_size(face,2)=stencil_size(face,2)+1;
            end
        end
        %
        % Adiciona as Faces ao Stencil %
        z=0;
        n=stencil_size(face,2);
        stencil_faces_aux=0;
        %
        for j=1:n
            %
            % Celula Selecionada %
            cell=stencil_cells(face,j);
            %
            % Verifica se a Celula Existe %
            if cell~=0
                for k=1:cell_face_num(cell)
                    %
                    % Seleciona uma Face da Celula %
                    face_aux=cell_faces(cell,k);
                    %
                    % Verifica se a Face é de Fronteira %
                    if face_bound(face_aux,1)==1
                        if z==0
                            z=z+1;
                            stencil_faces_aux(z)=face_aux;
                        else
                            teste=0;
                            for q=1:z
                                if stencil_faces_aux(q)==face_aux
                                    teste=1;
                                end
                            end
                            %
                            if teste==0
                                z=z+1;
                                stencil_faces_aux(z)=face_aux;
                            end
                        end
                    end
                end
            end
        end
        %
        % Organizar o Stencil das Faces %
        stencil_faces_aux=sort(stencil_faces_aux,'descend');
        stencil_size(face,3)=0;
        z=0;
        n=size(stencil_faces_aux,2);
        %
        for j=1:n
            if stencil_faces_aux~=0
                z=z+1;
                stencil_faces(face,z)=stencil_faces_aux(j);
                stencil_size(face,3)=stencil_size(face,3)+1;
            end
        end
        %
        % Tamanho do Stencil %
        stencil_size(face,1)=stencil_size(face,2)+stencil_size(face,3);
        %
        % Calcula o Tamanho do Stencil, Verifica se a Expansão Realizada foi Suficiente %
        x_max=0;
        x_min=L;
        y_max=0;
        y_min=L;
        %
        area_stencil=0;
        %
        % Células do Stencil %
        for j=1:stencil_size(face,2)
            %
            % Celula Selecionada %
            cell=stencil_cells(face,j);
            %
            % Comprimento em X %
            if x_max<cells(cell,1)
                x_max=cells(cell,1);
            end
            if x_min>cells(cell,1)
                x_min=cells(cell,1);
            end
            %
            % Comprimento em Y %
            if y_max<cells(cell,2)
                y_max=cells(cell,2);
            end
            if y_min>cells(cell,2)
                y_min=cells(cell,2);
            end
            %
            area_stencil=area_stencil+cell_vol(cell);
        end
        %
        % Faces do Stencil %
        for j=1:stencil_size(face,3)
            %
            % Face Selecionada %
            face_aux=stencil_faces(face,j);
            %
            % Comprimento em X %
            if x_max<faces(face_aux,1)
                x_max=faces(face_aux,1);
            end
            if x_min>faces(face_aux,1)
                x_min=faces(face_aux,1);
            end
            %
            % Comprimento em Y %
            if y_max<faces(face_aux,2)
                y_max=faces(face_aux,2);
            end
            if y_min>faces(face_aux,2)
                y_min=faces(face_aux,2);
            end
        end
        %
        % Determinação do Comprimento de Referencia do Stencil - Importante para os casos onde a malha não é
        % uniforme ou a malha é não estruturada %
        L_stencil=sqrt(area_stencil/stencil_size(face,2));
        %
        % Determina em que direcção o Stencil deve crescer %
        L_x=abs(x_max-x_min)/L_stencil;
        L_y=abs(y_max-y_min)/L_stencil;
        %
        % Determina quantos pontos tem o Stencil em cada direcção %
        n_x=int8(L_x+1.1);
        n_y=int8(L_y+1.1);
        %
    end
    stencil_size(face,1)=stencil_size(face,2)+stencil_size(face,3);
end
%incremento'
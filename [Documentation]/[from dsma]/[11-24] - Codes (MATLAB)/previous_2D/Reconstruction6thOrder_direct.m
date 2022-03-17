function [T,D]=Reconstruction6thOrder_direct(stencil_cells,stencil_faces,stencil_size,ponderado)
% warning('off','all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                  17 de Outubro de 2016                                            %
%                                  17 de Outubro de 2016                                            %
%                                                                                                   %
% Função que implementa a Solução Numérica a partir do Minimos Quadrados de 2ª Ordem                %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global L Lref cell_side vert_side cell_num face_num vert_num w;
global verts cells cell_verts cell_faces cell_vol faces face_area face_vert cell_norm;
global face_bound cell_bound face_cells vert_cells;
global phi lap_phi  phi_faces flux_phi_faces extended_stencil;

%
for i=1:face_num
    %
    face=i;
          % Normal Exterior da Face %
    vert=face_vert(face,:);
    nx=verts(vert(1),1)-verts(vert(2),1);
    ny=verts(vert(1),2)-verts(vert(2),2);
    
    %
    %% Construção da Matriz D para cada Face, Celulas do Stencil %%
    %
    n_cell=stencil_size(face,2);
    %
    for j=1:n_cell
        %
        % Celula Selecionada %
        cell=stencil_cells(face,j);
        
        % Matriz D %
        D{face}(j,1)=1;
        %
        D{face}(j,2)=(cells(cell,1)-faces(face,1));
        D{face}(j,3)=(cells(cell,2)-faces(face,2));
        %
        D{face}(j,4)=(cells(cell,1)-faces(face,1))^2;
        D{face}(j,5)=(cells(cell,2)-faces(face,2))^2;
        D{face}(j,6)=(cells(cell,1)-faces(face,1))*(cells(cell,2)-faces(face,2));
        %
        D{face}(j,7)=(cells(cell,1)-faces(face,1))^3;
        D{face}(j,8)=(cells(cell,2)-faces(face,2))^3;
        D{face}(j,9)=(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2));
        D{face}(j,10)=(cells(cell,1)-faces(face,1))*(cells(cell,2)-faces(face,2))^2;
        % 
        D{face}(j,11)=(cells(cell,1)-faces(face,1))^4;
        D{face}(j,12)=(cells(cell,2)-faces(face,2))^4;
        D{face}(j,13)=(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2));
        D{face}(j,14)=(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^2;
        D{face}(j,15)=(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^3;
        %
        D{face}(j,16)=(cells(cell,1)-faces(face,1))^5;
        D{face}(j,17)=(cells(cell,2)-faces(face,2))^5;
        D{face}(j,18)=(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^1;
        D{face}(j,19)=(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^2;
        D{face}(j,20)=(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^3;
        D{face}(j,21)=(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^4;
        %
        D{face}(j,22)=(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^1;
        D{face}(j,23)=(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^2;
        D{face}(j,24)=(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^3;
        D{face}(j,25)=(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^4;
        D{face}(j,26)=(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^5;
        %
        D{face}(j,27)=(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^2;
        D{face}(j,28)=(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^3;
        D{face}(j,29)=(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^4;
        D{face}(j,30)=(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^5;
        %
        D{face}(j,31)=(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^3;
        D{face}(j,32)=(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^4;
        D{face}(j,33)=(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^5;
        %
        D{face}(j,34)=(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^4;
        D{face}(j,35)=(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^5;
        %
        D{face}(j,36)=(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^5;
        
          
        if nx==0  
        D{face}(j,37)=(cells(cell,2)-faces(face,2))^6;
        D{face}(j,38)=(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^6;
        D{face}(j,39)=(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^6;
        D{face}(j,40)=(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^6;
        D{face}(j,41)=(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^6;
        D{face}(j,42)=(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^6;
        end
        
        if ny==0
        D{face}(j,37)=(cells(cell,1)-faces(face,1))^6;
        D{face}(j,38)=(cells(cell,1)-faces(face,1))^6*(cells(cell,2)-faces(face,2))^1;
        D{face}(j,39)=(cells(cell,1)-faces(face,1))^6*(cells(cell,2)-faces(face,2))^2;
        D{face}(j,40)=(cells(cell,1)-faces(face,1))^6*(cells(cell,2)-faces(face,2))^3;
        D{face}(j,41)=(cells(cell,1)-faces(face,1))^6*(cells(cell,2)-faces(face,2))^4;
        D{face}(j,42)=(cells(cell,1)-faces(face,1))^6*(cells(cell,2)-faces(face,2))^5;
        end
        
    end
    %
    %% Construção da Matriz D para cada Face, Faces do Stencil %%
    %
    n_face=stencil_size(face,3);
    %
    for j=1:n_face
        %
        % Face Selecionada %
        facestencil=stencil_faces(face,j);
        %
        
        % Matriz D %
        D{face}(n_cell+j,1)=1;
        %
        D{face}(n_cell+j,2)=(faces(facestencil,1)-faces(face,1));
        D{face}(n_cell+j,3)=(faces(facestencil,2)-faces(face,2));
        %
        D{face}(n_cell+j,4)=(faces(facestencil,1)-faces(face,1))^2;
        D{face}(n_cell+j,5)=(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,6)=(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2));
        %
        D{face}(n_cell+j,7)=(faces(facestencil,1)-faces(face,1))^3;
        D{face}(n_cell+j,8)=(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,9)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2));
        D{face}(n_cell+j,10)=(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2;
        %
        D{face}(n_cell+j,11)=(faces(facestencil,1)-faces(face,1))^4;
        D{face}(n_cell+j,12)=(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,13)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,14)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,15)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3;
        %
        D{face}(n_cell+j,16)=(faces(facestencil,1)-faces(face,1))^5;
        D{face}(n_cell+j,17)=(faces(facestencil,2)-faces(face,2))^5;
        D{face}(n_cell+j,18)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,19)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,20)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,21)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4;
        %
        D{face}(n_cell+j,22)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,23)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,24)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,25)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,26)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5;
        %
        D{face}(n_cell+j,27)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,28)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,29)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,30)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^5;
        %
        D{face}(n_cell+j,31)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,32)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,33)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^5;
        %
        D{face}(n_cell+j,34)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,35)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^5;
        %
        D{face}(n_cell+j,36)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^5;
        %
        
        if nx==0  
        D{face}(n_cell+j,37)=(faces(facestencil,2)-faces(face,2))^6;
        D{face}(n_cell+j,38)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^6;
        D{face}(n_cell+j,39)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^6;
        D{face}(n_cell+j,40)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^6;
        D{face}(n_cell+j,41)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^6;
        D{face}(n_cell+j,42)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^6;
        end
        
        if ny==0
        D{face}(n_cell+j,37)=(faces(facestencil,1)-faces(face,1))^6;
        D{face}(n_cell+j,38)=(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,39)=(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,40)=(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,41)=(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,42)=(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^5;
        end
        
        
        
    end
    %
end
    %% Determinação da Matriz T %%
for face=1:face_num
    
    % D*c=phi -> c=inv(D'D)*D'*phi -> T=inv(D'D)*D' %
      T{face}=inv(D{face});
%     T{face}=matrix_inverter(D{face}'*D1{face})*D1{face}';
    
end
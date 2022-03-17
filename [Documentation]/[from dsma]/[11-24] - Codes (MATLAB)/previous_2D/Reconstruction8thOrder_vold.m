function [T,D]=Reconstruction8thOrder(stencil_cells,stencil_faces,stencil_size,ponderado)
% warning('off','all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                  17 de Outubro de 2016                                            %
%                                  17 de Outubro de 2016                                            %
%                                                                                                   %
% Fun��o que implementa a Solu��o Num�rica a partir do Minimos Quadrados de 2� Ordem                %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global L Lref cell_side vert_side cell_num face_num vert_num w;
global verts cells cell_verts cell_faces cell_vol faces face_area face_vert cell_norm;
global face_bound cell_bound face_cells vert_cells u_convec_x u_convec_y gamma_diff;
global phi lap_phi  phi_faces flux_phi_faces extended_stencil dimensional_correction neuman robin;

%
for i=1:face_num
    %
    face=i;
          % Normal Exterior da Face %
    vert=face_vert(face,:);
    nx=verts(vert(1),1)-verts(vert(2),1);
    ny=verts(vert(1),2)-verts(vert(2),2);
    
    %
    %% Constru��o da Matriz D para cada Face, Celulas do Stencil %%
    %
    n_cell=stencil_size(face,2);
    %
    for j=1:n_cell
        %
        % Celula Selecionada %
        cell=stencil_cells(face,j);
        %
        if ponderado==true
            dx=(cells(cell,1)-faces(face,1));
            dy=(cells(cell,2)-faces(face,2));
            w=1/(sqrt(dx^2+dy^2))^1;
        else
            w=1;
        end
        %
        % Matriz D tendo em conta a pondera��o %
        D1{face}(j,1)=w*1;
        %
        D1{face}(j,2)=w*(cells(cell,1)-faces(face,1));
        D1{face}(j,3)=w*(cells(cell,2)-faces(face,2));
        %
        D1{face}(j,4)=w*(cells(cell,1)-faces(face,1))^2;
        D1{face}(j,5)=w*(cells(cell,2)-faces(face,2))^2;
        D1{face}(j,6)=w*(cells(cell,1)-faces(face,1))*(cells(cell,2)-faces(face,2));
        %
        D1{face}(j,7)=w*(cells(cell,1)-faces(face,1))^3;
        D1{face}(j,8)=w*(cells(cell,2)-faces(face,2))^3;
        D1{face}(j,9)=w*(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2));
        D1{face}(j,10)=w*(cells(cell,1)-faces(face,1))*(cells(cell,2)-faces(face,2))^2;
        % 
        D1{face}(j,11)=w*(cells(cell,1)-faces(face,1))^4;
        D1{face}(j,12)=w*(cells(cell,2)-faces(face,2))^4;
        D1{face}(j,13)=w*(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2));
        D1{face}(j,14)=w*(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^2;
        D1{face}(j,15)=w*(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^3;
        %
        D1{face}(j,16)=w*(cells(cell,1)-faces(face,1))^5;
        D1{face}(j,17)=w*(cells(cell,2)-faces(face,2))^5;
        D1{face}(j,18)=w*(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^1;
        D1{face}(j,19)=w*(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^2;
        D1{face}(j,20)=w*(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^3;
        D1{face}(j,21)=w*(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^4;
        %
        D1{face}(j,22)=w*(cells(cell,1)-faces(face,1))^6;
        D1{face}(j,23)=w*(cells(cell,2)-faces(face,2))^6;
        D1{face}(j,24)=w*(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^1;
        D1{face}(j,25)=w*(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^2;
        D1{face}(j,26)=w*(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^3;
        D1{face}(j,27)=w*(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^4;
        D1{face}(j,28)=w*(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^5;
        %
        D1{face}(j,29)=w*(cells(cell,1)-faces(face,1))^7;
        D1{face}(j,30)=w*(cells(cell,2)-faces(face,2))^7;
        D1{face}(j,31)=w*(cells(cell,1)-faces(face,1))^6*(cells(cell,2)-faces(face,2))^1;
        D1{face}(j,32)=w*(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^2;
        D1{face}(j,33)=w*(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^3;
        D1{face}(j,34)=w*(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^4;
        D1{face}(j,35)=w*(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^5;
        D1{face}(j,36)=w*(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^6;
        %
        
        if extended_stencil && nx==0       
        D1{face}(j,37)=w*(cells(cell,2)-faces(face,2))^8;
        D1{face}(j,38)=w*(cells(cell,1)-faces(face,1))^7*(cells(cell,2)-faces(face,2))^1;
        D1{face}(j,39)=w*(cells(cell,1)-faces(face,1))^6*(cells(cell,2)-faces(face,2))^2;
        D1{face}(j,40)=w*(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^3;
        D1{face}(j,41)=w*(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^4;
        D1{face}(j,42)=w*(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^5;
        D1{face}(j,43)=w*(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^6;
        D1{face}(j,44)=w*(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^7;
        elseif extended_stencil && ny==0
        D1{face}(j,37)=w*(cells(cell,1)-faces(face,1))^8;
        D1{face}(j,38)=w*(cells(cell,1)-faces(face,1))^7*(cells(cell,2)-faces(face,2))^1;
        D1{face}(j,39)=w*(cells(cell,1)-faces(face,1))^6*(cells(cell,2)-faces(face,2))^2;
        D1{face}(j,40)=w*(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^3;
        D1{face}(j,41)=w*(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^4;
        D1{face}(j,42)=w*(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^5;
        D1{face}(j,43)=w*(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^6;
        D1{face}(j,44)=w*(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^7;
        end
        
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
        D{face}(j,22)=(cells(cell,1)-faces(face,1))^6;
        D{face}(j,23)=(cells(cell,2)-faces(face,2))^6;
        D{face}(j,24)=(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^1;
        D{face}(j,25)=(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^2;
        D{face}(j,26)=(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^3;
        D{face}(j,27)=(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^4;
        D{face}(j,28)=(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^5;
        %
        D{face}(j,29)=(cells(cell,1)-faces(face,1))^7;
        D{face}(j,30)=(cells(cell,2)-faces(face,2))^7;
        D{face}(j,31)=(cells(cell,1)-faces(face,1))^6*(cells(cell,2)-faces(face,2))^1;
        D{face}(j,32)=(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^2;
        D{face}(j,33)=(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^3;
        D{face}(j,34)=(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^4;
        D{face}(j,35)=(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^5;
        D{face}(j,36)=(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^6;
        
        
        if extended_stencil && nx==0       
        D{face}(j,37)=(cells(cell,2)-faces(face,2))^8;
        D{face}(j,38)=(cells(cell,1)-faces(face,1))^7*(cells(cell,2)-faces(face,2))^1;
        D{face}(j,39)=(cells(cell,1)-faces(face,1))^6*(cells(cell,2)-faces(face,2))^2;
        D{face}(j,40)=(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^3;
        D{face}(j,41)=(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^4;
        D{face}(j,42)=(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^5;
        D{face}(j,43)=(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^6;
        D{face}(j,44)=(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^7;
        elseif extended_stencil && ny==0
        D{face}(j,37)=(cells(cell,1)-faces(face,1))^8;
        D{face}(j,38)=(cells(cell,1)-faces(face,1))^7*(cells(cell,2)-faces(face,2))^1;
        D{face}(j,39)=(cells(cell,1)-faces(face,1))^6*(cells(cell,2)-faces(face,2))^2;
        D{face}(j,40)=(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^3;
        D{face}(j,41)=(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^4;
        D{face}(j,42)=(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^5;
        D{face}(j,43)=(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^6;
        D{face}(j,44)=(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^7;
        end
        
    end
    %
    %% Constru��o da Matriz D para cada Face, Faces do Stencil %%
    %
    n_face=stencil_size(face,3);
    %
    for j=1:n_face
        %
        % Face Selecionada %
        facestencil=stencil_faces(face,j);
        %
        if ponderado==true
            if facestencil~=face
                dx=(faces(facestencil,1)-faces(face,1));
                dy=(faces(facestencil,2)-faces(face,2));
            else
                aux=face_cells(face,1);
                dx=(cells(aux,1)-faces(face,1));
                dy=(cells(aux,2)-faces(face,2));
%                      dx=0.0000001;
%                     dy=0.0000001;

            end
            w=1/(sqrt(dx^2+dy^2))^1;
        else
            w=1;
        end
        %
        % Matriz D tendo em conta a pondera��o %
        %% NEUMAN     
 
 
            normal_x=face_bound(facestencil,2)*face_area(facestencil,1);
            normal_y=face_bound(facestencil,3)*face_area(facestencil,1);
            
%             if facestencil<cell_side+1
%             normal_x=-1*face_area(facestencil,1);    
%             end
            
            if ~dimensional_correction
                      normal_x=normal_x/face_area(facestencil,1);
                      normal_y=normal_y/face_area(facestencil,1);
              end 
            
            
            
            if neuman && face_bound(facestencil,4)==1
                
        D1{face}(n_cell+j,1)=0;
        %
        D1{face}(n_cell+j,2)=w*normal_x;
        D1{face}(n_cell+j,3)=0;
        %
        D1{face}(n_cell+j,4)=2*w*(faces(facestencil,1)-faces(face,1))*normal_x;
        D1{face}(n_cell+j,5)=0;
        D1{face}(n_cell+j,6)=w*(faces(facestencil,2)-faces(face,2))*normal_x;
        %
        D1{face}(n_cell+j,7)=3*w*(faces(facestencil,1)-faces(face,1))^2*normal_x;
        D1{face}(n_cell+j,8)=0;
        D1{face}(n_cell+j,9)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))*normal_x;
        D1{face}(n_cell+j,10)=w*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        %
        D1{face}(n_cell+j,11)=4*w*(faces(facestencil,1)-faces(face,1))^3*normal_x;
        D1{face}(n_cell+j,12)=0;
        D1{face}(n_cell+j,13)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,14)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,15)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        %
        D1{face}(n_cell+j,16)=5*w*(faces(facestencil,1)-faces(face,1))^4*normal_x;
        D1{face}(n_cell+j,17)=0;
        D1{face}(n_cell+j,18)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,19)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,20)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D1{face}(n_cell+j,21)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        %
        D1{face}(n_cell+j,22)=6*w*(faces(facestencil,1)-faces(face,1))^5*normal_x;
        D1{face}(n_cell+j,23)=0;
        D1{face}(n_cell+j,24)=5*w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,25)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,26)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D1{face}(n_cell+j,27)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D1{face}(n_cell+j,28)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        %
        D1{face}(n_cell+j,29)=7*w*(faces(facestencil,1)-faces(face,1))^6*normal_x;
        D1{face}(n_cell+j,30)=0;
        D1{face}(n_cell+j,31)=6*w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,32)=5*w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,33)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D1{face}(n_cell+j,34)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D1{face}(n_cell+j,35)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        D1{face}(n_cell+j,36)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^6*normal_x;
        %
        
        if extended_stencil && nx==0       
        D1{face}(n_cell+j,37)=0;
        D1{face}(n_cell+j,38)=7*w*(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,39)=6*w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,40)=5*w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D1{face}(n_cell+j,41)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D1{face}(n_cell+j,42)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        D1{face}(n_cell+j,43)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^6*normal_x;
        D1{face}(n_cell+j,44)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^7*normal_x;
        elseif extended_stencil && ny==0
        D1{face}(n_cell+j,37)=8*w*(faces(facestencil,1)-faces(face,1))^7*normal_x;
        D1{face}(n_cell+j,38)=7*w*(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,39)=6*w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,40)=5*w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D1{face}(n_cell+j,41)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D1{face}(n_cell+j,42)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        D1{face}(n_cell+j,43)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^6*normal_x;
        D1{face}(n_cell+j,44)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^7*normal_x;
        end
        
        % Matriz D %
        D{face}(n_cell+j,1)=0;
        %
        D{face}(n_cell+j,2)=normal_x;
        D{face}(n_cell+j,3)=0;
        %
        D{face}(n_cell+j,4)=2*(faces(facestencil,1)-faces(face,1))*normal_x;
        D{face}(n_cell+j,5)=0;
        D{face}(n_cell+j,6)=(faces(facestencil,2)-faces(face,2))*normal_x;
        %
        D{face}(n_cell+j,7)=3*(faces(facestencil,1)-faces(face,1))^2*normal_x;
        D{face}(n_cell+j,8)=0;
        D{face}(n_cell+j,9)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))*normal_x;
        D{face}(n_cell+j,10)=(faces(facestencil,2)-faces(face,2))^2*normal_x;
        %
        D{face}(n_cell+j,11)=4*(faces(facestencil,1)-faces(face,1))^3*normal_x;
        D{face}(n_cell+j,12)=0;
        D{face}(n_cell+j,13)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,14)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,15)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        %
        D{face}(n_cell+j,16)=5*(faces(facestencil,1)-faces(face,1))^4*normal_x;
        D{face}(n_cell+j,17)=0;
        D{face}(n_cell+j,18)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,19)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,20)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D{face}(n_cell+j,21)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        %
        D{face}(n_cell+j,22)=6*(faces(facestencil,1)-faces(face,1))^5*normal_x;
        D{face}(n_cell+j,23)=0;
        D{face}(n_cell+j,24)=5*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,25)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,26)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D{face}(n_cell+j,27)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D{face}(n_cell+j,28)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        %
        D{face}(n_cell+j,29)=7*(faces(facestencil,1)-faces(face,1))^6*normal_x;
        D{face}(n_cell+j,30)=0;
        D{face}(n_cell+j,31)=6*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,32)=5*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,33)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D{face}(n_cell+j,34)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D{face}(n_cell+j,35)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        D{face}(n_cell+j,36)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^6*normal_x;
        %
        
        if extended_stencil && nx==0       
        D{face}(n_cell+j,37)=0;
        D{face}(n_cell+j,38)=7*(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,39)=6*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,40)=5*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D{face}(n_cell+j,41)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D{face}(n_cell+j,42)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        D{face}(n_cell+j,43)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^6*normal_x;
        D{face}(n_cell+j,44)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^7*normal_x;
        elseif extended_stencil && ny==0
        D{face}(n_cell+j,37)=8*(faces(facestencil,1)-faces(face,1))^7*normal_x;
        D{face}(n_cell+j,38)=7*(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,39)=6*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,40)=5*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D{face}(n_cell+j,41)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D{face}(n_cell+j,42)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        D{face}(n_cell+j,43)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^6*normal_x;
        D{face}(n_cell+j,44)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^7*normal_x;
        end
%% ROBIN



        elseif robin && face_bound(facestencil,4)==2
                
        D1{face}(n_cell+j,1)=0;
        %
        D1{face}(n_cell+j,2)=w*normal_x;
        D1{face}(n_cell+j,3)=0;
        %
        D1{face}(n_cell+j,4)=2*w*(faces(facestencil,1)-faces(face,1))*normal_x;
        D1{face}(n_cell+j,5)=0;
        D1{face}(n_cell+j,6)=w*(faces(facestencil,2)-faces(face,2))*normal_x;
        %
        D1{face}(n_cell+j,7)=3*w*(faces(facestencil,1)-faces(face,1))^2*normal_x;
        D1{face}(n_cell+j,8)=0;
        D1{face}(n_cell+j,9)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))*normal_x;
        D1{face}(n_cell+j,10)=w*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        %
        D1{face}(n_cell+j,11)=4*w*(faces(facestencil,1)-faces(face,1))^3*normal_x;
        D1{face}(n_cell+j,12)=0;
        D1{face}(n_cell+j,13)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,14)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,15)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        %
        D1{face}(n_cell+j,16)=5*w*(faces(facestencil,1)-faces(face,1))^4*normal_x;
        D1{face}(n_cell+j,17)=0;
        D1{face}(n_cell+j,18)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,19)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,20)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D1{face}(n_cell+j,21)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        %
        D1{face}(n_cell+j,22)=6*w*(faces(facestencil,1)-faces(face,1))^5*normal_x;
        D1{face}(n_cell+j,23)=0;
        D1{face}(n_cell+j,24)=5*w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,25)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,26)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D1{face}(n_cell+j,27)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D1{face}(n_cell+j,28)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        %
        D1{face}(n_cell+j,29)=7*w*(faces(facestencil,1)-faces(face,1))^6*normal_x;
        D1{face}(n_cell+j,30)=0;
        D1{face}(n_cell+j,31)=6*w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,32)=5*w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,33)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D1{face}(n_cell+j,34)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D1{face}(n_cell+j,35)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        D1{face}(n_cell+j,36)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^6*normal_x;
        %
        
        if extended_stencil && nx==0       
        D1{face}(n_cell+j,37)=0;
        D1{face}(n_cell+j,38)=7*w*(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,39)=6*w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,40)=5*w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D1{face}(n_cell+j,41)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D1{face}(n_cell+j,42)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        D1{face}(n_cell+j,43)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^6*normal_x;
        D1{face}(n_cell+j,44)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^7*normal_x;
        elseif extended_stencil && ny==0
        D1{face}(n_cell+j,37)=8*w*(faces(facestencil,1)-faces(face,1))^7*normal_x;
        D1{face}(n_cell+j,38)=7*w*(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,39)=6*w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,40)=5*w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D1{face}(n_cell+j,41)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D1{face}(n_cell+j,42)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        D1{face}(n_cell+j,43)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^6*normal_x;
        D1{face}(n_cell+j,44)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^7*normal_x;
        end
        
        % Matriz D %
        D{face}(n_cell+j,1)=0;
        %
        D{face}(n_cell+j,2)=normal_x;
        D{face}(n_cell+j,3)=0;
        %
        D{face}(n_cell+j,4)=2*(faces(facestencil,1)-faces(face,1))*normal_x;
        D{face}(n_cell+j,5)=0;
        D{face}(n_cell+j,6)=(faces(facestencil,2)-faces(face,2))*normal_x;
        %
        D{face}(n_cell+j,7)=3*(faces(facestencil,1)-faces(face,1))^2*normal_x;
        D{face}(n_cell+j,8)=0;
        D{face}(n_cell+j,9)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))*normal_x;
        D{face}(n_cell+j,10)=(faces(facestencil,2)-faces(face,2))^2*normal_x;
        %
        D{face}(n_cell+j,11)=4*(faces(facestencil,1)-faces(face,1))^3*normal_x;
        D{face}(n_cell+j,12)=0;
        D{face}(n_cell+j,13)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,14)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,15)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        %
        D{face}(n_cell+j,16)=5*(faces(facestencil,1)-faces(face,1))^4*normal_x;
        D{face}(n_cell+j,17)=0;
        D{face}(n_cell+j,18)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,19)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,20)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D{face}(n_cell+j,21)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        %
        D{face}(n_cell+j,22)=6*(faces(facestencil,1)-faces(face,1))^5*normal_x;
        D{face}(n_cell+j,23)=0;
        D{face}(n_cell+j,24)=5*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,25)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,26)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D{face}(n_cell+j,27)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D{face}(n_cell+j,28)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        %
        D{face}(n_cell+j,29)=7*(faces(facestencil,1)-faces(face,1))^6*normal_x;
        D{face}(n_cell+j,30)=0;
        D{face}(n_cell+j,31)=6*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,32)=5*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,33)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D{face}(n_cell+j,34)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D{face}(n_cell+j,35)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        D{face}(n_cell+j,36)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^6*normal_x;
        %
        
        if extended_stencil && nx==0       
        D{face}(n_cell+j,37)=0;
        D{face}(n_cell+j,38)=7*(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,39)=6*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,40)=5*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D{face}(n_cell+j,41)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D{face}(n_cell+j,42)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        D{face}(n_cell+j,43)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^6*normal_x;
        D{face}(n_cell+j,44)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^7*normal_x;
        elseif extended_stencil && ny==0
        D{face}(n_cell+j,37)=8*(faces(facestencil,1)-faces(face,1))^7*normal_x;
        D{face}(n_cell+j,38)=7*(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,39)=6*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,40)=5*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D{face}(n_cell+j,41)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D{face}(n_cell+j,42)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        D{face}(n_cell+j,43)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^6*normal_x;
        D{face}(n_cell+j,44)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^7*normal_x;
        end
%%%%%%%%%

        D1{face}(n_cell+j,1)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,1)+w*1;
        %
        D1{face}(n_cell+j,2)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,2)+w*(faces(facestencil,1)-faces(face,1));
        D1{face}(n_cell+j,3)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,3)+w*(faces(facestencil,2)-faces(face,2));
        %
        D1{face}(n_cell+j,4)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,4)+w*(faces(facestencil,1)-faces(face,1))^2;
        D1{face}(n_cell+j,5)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,5)+w*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,6)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,6)+w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2));
        %
        D1{face}(n_cell+j,7)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,7)+w*(faces(facestencil,1)-faces(face,1))^3;
        D1{face}(n_cell+j,8)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,8)+w*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,9)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,9)+w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2));
        D1{face}(n_cell+j,10)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,10)+w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2;
        %
        D1{face}(n_cell+j,11)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,11)+w*(faces(facestencil,1)-faces(face,1))^4;
        D1{face}(n_cell+j,12)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,12)+w*(faces(facestencil,2)-faces(face,2))^4;
        D1{face}(n_cell+j,13)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,13)+w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,14)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,14)+w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,15)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,15)+w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3;
        %
        D1{face}(n_cell+j,16)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,16)+w*(faces(facestencil,1)-faces(face,1))^5;
        D1{face}(n_cell+j,17)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,17)+w*(faces(facestencil,2)-faces(face,2))^5;
        D1{face}(n_cell+j,18)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,18)+w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,19)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,19)+w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,20)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,20)+w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,21)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,21)+w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4;
        %
        D1{face}(n_cell+j,22)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,22)+w*(faces(facestencil,1)-faces(face,1))^6;
        D1{face}(n_cell+j,23)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,23)+w*(faces(facestencil,2)-faces(face,2))^6;
        D1{face}(n_cell+j,24)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,24)+w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,25)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,25)+w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,26)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,26)+w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,27)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,27)+w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4;
        D1{face}(n_cell+j,28)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,28)+w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5;
        %
        D1{face}(n_cell+j,29)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,29)+w*(faces(facestencil,1)-faces(face,1))^7;
        D1{face}(n_cell+j,30)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,30)+w*(faces(facestencil,2)-faces(face,2))^7;
        D1{face}(n_cell+j,31)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,31)+w*(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,32)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,32)+w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,33)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,33)+w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,34)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,34)+w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^4;
        D1{face}(n_cell+j,35)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,35)+w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^5;
        D1{face}(n_cell+j,36)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,36)+w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^6;
        %
        
        if extended_stencil && nx==0       
        D1{face}(n_cell+j,37)=w*(faces(facestencil,2)-faces(face,2))^8;
        D1{face}(n_cell+j,38)=w*(faces(facestencil,1)-faces(face,1))^7*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,39)=w*(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,40)=w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,41)=w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^4;
        D1{face}(n_cell+j,42)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^5;
        D1{face}(n_cell+j,43)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^6;
        D1{face}(n_cell+j,44)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^7;
        elseif extended_stencil && ny==0
        D1{face}(n_cell+j,37)=w*(faces(facestencil,1)-faces(face,1))^8;
        D1{face}(n_cell+j,38)=w*(faces(facestencil,1)-faces(face,1))^7*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,39)=w*(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,40)=w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,41)=w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^4;
        D1{face}(n_cell+j,42)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^5;
        D1{face}(n_cell+j,43)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^6;
        D1{face}(n_cell+j,44)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^7;
        end
        
        % Matriz D %
        D{face}(n_cell+j,1)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,1)+1;
        %
        D{face}(n_cell+j,2)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,2)+(faces(facestencil,1)-faces(face,1));
        D{face}(n_cell+j,3)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,3)+(faces(facestencil,2)-faces(face,2));
        %
        D{face}(n_cell+j,4)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,4)+(faces(facestencil,1)-faces(face,1))^2;
        D{face}(n_cell+j,5)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,5)+(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,6)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,6)+(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2));
        %
        D{face}(n_cell+j,7)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,7)+(faces(facestencil,1)-faces(face,1))^3;
        D{face}(n_cell+j,8)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,8)+(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,9)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,9)+(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2));
        D{face}(n_cell+j,10)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,10)+(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2;
        %
        D{face}(n_cell+j,11)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,11)+(faces(facestencil,1)-faces(face,1))^4;
        D{face}(n_cell+j,12)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,12)+(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,13)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,13)+(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,14)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,14)+(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,15)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,15)+(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3;
        %
        D{face}(n_cell+j,16)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,16)+(faces(facestencil,1)-faces(face,1))^5;
        D{face}(n_cell+j,17)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,17)+(faces(facestencil,2)-faces(face,2))^5;
        D{face}(n_cell+j,18)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,18)+(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,19)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,19)+(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,20)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,20)+(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,21)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,21)+(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4;
        %
        D{face}(n_cell+j,22)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,22)+(faces(facestencil,1)-faces(face,1))^6;
        D{face}(n_cell+j,23)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,23)+(faces(facestencil,2)-faces(face,2))^6;
        D{face}(n_cell+j,24)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,24)+(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,25)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,25)+(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,26)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,26)+(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,27)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,27)+(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,28)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,28)+(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5;
        %
        D{face}(n_cell+j,29)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,29)+(faces(facestencil,1)-faces(face,1))^7;
        D{face}(n_cell+j,30)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,30)+(faces(facestencil,2)-faces(face,2))^7;
        D{face}(n_cell+j,31)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,31)+(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,32)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,32)+(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,33)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,33)+(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,34)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,34)+(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,35)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,35)+(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^5;
        D{face}(n_cell+j,36)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,36)+(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^6;
        
        
        
        if extended_stencil && nx==0       
        D{face}(n_cell+j,37)=(faces(facestencil,2)-faces(face,2))^8;
        D{face}(n_cell+j,38)=(faces(facestencil,1)-faces(face,1))^7*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,39)=(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,40)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,41)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,42)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^5;
        D{face}(n_cell+j,43)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^6;
        D{face}(n_cell+j,44)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^7;
        elseif extended_stencil && ny==0
        D{face}(n_cell+j,37)=(faces(facestencil,1)-faces(face,1))^8;
        D{face}(n_cell+j,38)=(faces(facestencil,1)-faces(face,1))^7*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,39)=(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,40)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,41)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,42)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^5;
        D{face}(n_cell+j,43)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^6;
        D{face}(n_cell+j,44)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^7;
        end



                
%% Dirichlet                
            else
        D1{face}(n_cell+j,1)=w*1;
        %
        D1{face}(n_cell+j,2)=w*(faces(facestencil,1)-faces(face,1));
        D1{face}(n_cell+j,3)=w*(faces(facestencil,2)-faces(face,2));
        %
        D1{face}(n_cell+j,4)=w*(faces(facestencil,1)-faces(face,1))^2;
        D1{face}(n_cell+j,5)=w*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,6)=w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2));
        %
        D1{face}(n_cell+j,7)=w*(faces(facestencil,1)-faces(face,1))^3;
        D1{face}(n_cell+j,8)=w*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,9)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2));
        D1{face}(n_cell+j,10)=w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2;
        %
        D1{face}(n_cell+j,11)=w*(faces(facestencil,1)-faces(face,1))^4;
        D1{face}(n_cell+j,12)=w*(faces(facestencil,2)-faces(face,2))^4;
        D1{face}(n_cell+j,13)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,14)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,15)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3;
        %
        D1{face}(n_cell+j,16)=w*(faces(facestencil,1)-faces(face,1))^5;
        D1{face}(n_cell+j,17)=w*(faces(facestencil,2)-faces(face,2))^5;
        D1{face}(n_cell+j,18)=w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,19)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,20)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,21)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4;
        %
        D1{face}(n_cell+j,22)=w*(faces(facestencil,1)-faces(face,1))^6;
        D1{face}(n_cell+j,23)=w*(faces(facestencil,2)-faces(face,2))^6;
        D1{face}(n_cell+j,24)=w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,25)=w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,26)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,27)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4;
        D1{face}(n_cell+j,28)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5;
        %
        D1{face}(n_cell+j,29)=w*(faces(facestencil,1)-faces(face,1))^7;
        D1{face}(n_cell+j,30)=w*(faces(facestencil,2)-faces(face,2))^7;
        D1{face}(n_cell+j,31)=w*(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,32)=w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,33)=w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,34)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^4;
        D1{face}(n_cell+j,35)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^5;
        D1{face}(n_cell+j,36)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^6;
        %
        
        if extended_stencil && nx==0       
        D1{face}(n_cell+j,37)=w*(faces(facestencil,2)-faces(face,2))^8;
        D1{face}(n_cell+j,38)=w*(faces(facestencil,1)-faces(face,1))^7*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,39)=w*(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,40)=w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,41)=w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^4;
        D1{face}(n_cell+j,42)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^5;
        D1{face}(n_cell+j,43)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^6;
        D1{face}(n_cell+j,44)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^7;
        elseif extended_stencil && ny==0
        D1{face}(n_cell+j,37)=w*(faces(facestencil,1)-faces(face,1))^8;
        D1{face}(n_cell+j,38)=w*(faces(facestencil,1)-faces(face,1))^7*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,39)=w*(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,40)=w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,41)=w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^4;
        D1{face}(n_cell+j,42)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^5;
        D1{face}(n_cell+j,43)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^6;
        D1{face}(n_cell+j,44)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^7;
        end
        
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
        D{face}(n_cell+j,22)=(faces(facestencil,1)-faces(face,1))^6;
        D{face}(n_cell+j,23)=(faces(facestencil,2)-faces(face,2))^6;
        D{face}(n_cell+j,24)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,25)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,26)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,27)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,28)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5;
        %
        D{face}(n_cell+j,29)=(faces(facestencil,1)-faces(face,1))^7;
        D{face}(n_cell+j,30)=(faces(facestencil,2)-faces(face,2))^7;
        D{face}(n_cell+j,31)=(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,32)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,33)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,34)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,35)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^5;
        D{face}(n_cell+j,36)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^6;
        
        
        
        if extended_stencil && nx==0       
        D{face}(n_cell+j,37)=(faces(facestencil,2)-faces(face,2))^8;
        D{face}(n_cell+j,38)=(faces(facestencil,1)-faces(face,1))^7*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,39)=(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,40)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,41)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,42)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^5;
        D{face}(n_cell+j,43)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^6;
        D{face}(n_cell+j,44)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^7;
        elseif extended_stencil && ny==0
        D{face}(n_cell+j,37)=(faces(facestencil,1)-faces(face,1))^8;
        D{face}(n_cell+j,38)=(faces(facestencil,1)-faces(face,1))^7*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,39)=(faces(facestencil,1)-faces(face,1))^6*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,40)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,41)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,42)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^5;
        D{face}(n_cell+j,43)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^6;
        D{face}(n_cell+j,44)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^7;
        end
            end
        
    end
    %
end
    %% Determina��o da Matriz T %%
for face=1:face_num
    
    % D*c=phi -> c=inv(D'D)*D'*phi -> T=inv(D'D)*D' %
      T{face}=inv(D{face}'*D1{face})*D1{face}';
%     T{face}=matrix_inverter(D{face}'*D1{face})*D1{face}';
    
end
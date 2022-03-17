function [T,D]=Reconstruction6thOrder(stencil_cells,stencil_faces,stencil_size,ponderado)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                  13 de Outubro de 2016                                            %
%                                  14 de Outubro de 2016                                            %
%                                                                                                   %
% Função que implementa a Solução Numérica a partir do Minimos Quadrados de 2ª Ordem                %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global L Lref cell_side vert_side cell_num face_num vert_num w;
global verts cells cell_verts cell_faces cell_vol faces face_area face_vert cell_norm;
global face_bound cell_bound face_cells vert_cells u_convec_x u_convec_y gamma_diff;
global phi lap_phi  phi_faces flux_phi_faces face_w_gauss extended_stencil neuman dimensional_correction robin;


inteiro=0;
for i=1:face_num
    %
    face=i;
    
      % Normal Exterior da Face %
    vert=face_vert(face,:);
    nx=verts(vert(1),1)-verts(vert(2),1);
    ny=verts(vert(1),2)-verts(vert(2),2);
            
    %
    %% Construção da Matriz D para cada Face, Celulas do Stencil %%
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
            w=1/(sqrt(dx^2+dy^2))^6;
        else
            w=1;
        end
        %
        % Matriz D tendo em conta a ponderação %
        D1{face}(j,1)=w*1;
        %
        D1{face}(j,2)=w*(cells(cell,1)-faces(face,1)+inteiro);
        D1{face}(j,3)=w*(cells(cell,2)-faces(face,2)+inteiro);
        %
        D1{face}(j,4)=w*(cells(cell,1)-faces(face,1)+inteiro)^2;
        D1{face}(j,5)=w*(cells(cell,2)-faces(face,2)+inteiro)^2;
        D1{face}(j,6)=w*(cells(cell,1)-faces(face,1)+inteiro)*(cells(cell,2)-faces(face,2)+inteiro);
        %
        D1{face}(j,7)=w*(cells(cell,1)-faces(face,1)+inteiro)^3;
        D1{face}(j,8)=w*(cells(cell,2)-faces(face,2)+inteiro)^3;
        D1{face}(j,9)=w*(cells(cell,1)-faces(face,1)+inteiro)^2*(cells(cell,2)-faces(face,2)+inteiro);
        D1{face}(j,10)=w*(cells(cell,1)-faces(face,1)+inteiro)*(cells(cell,2)-faces(face,2)+inteiro)^2;
        % 
        D1{face}(j,11)=w*(cells(cell,1)-faces(face,1)+inteiro)^4;
        D1{face}(j,12)=w*(cells(cell,2)-faces(face,2)+inteiro)^4;
        D1{face}(j,13)=w*(cells(cell,1)-faces(face,1)+inteiro)^3*(cells(cell,2)-faces(face,2)+inteiro);
        D1{face}(j,14)=w*(cells(cell,1)-faces(face,1)+inteiro)^2*(cells(cell,2)-faces(face,2)+inteiro)^2;
        D1{face}(j,15)=w*(cells(cell,1)-faces(face,1)+inteiro)^1*(cells(cell,2)-faces(face,2)+inteiro)^3;
        %
        D1{face}(j,16)=w*(cells(cell,1)-faces(face,1)+inteiro)^5;
        D1{face}(j,17)=w*(cells(cell,2)-faces(face,2)+inteiro)^5;
        D1{face}(j,18)=w*(cells(cell,1)-faces(face,1)+inteiro)^4*(cells(cell,2)-faces(face,2)+inteiro)^1;
        D1{face}(j,19)=w*(cells(cell,1)-faces(face,1)+inteiro)^3*(cells(cell,2)-faces(face,2)+inteiro)^2;
        D1{face}(j,20)=w*(cells(cell,1)-faces(face,1)+inteiro)^2*(cells(cell,2)-faces(face,2)+inteiro)^3;
        D1{face}(j,21)=w*(cells(cell,1)-faces(face,1)+inteiro)^1*(cells(cell,2)-faces(face,2)+inteiro)^4;
        %
        
        if extended_stencil && nx==0       
        D1{face}(j,22)=w*(cells(cell,2)-faces(face,2))^6;
        D1{face}(j,23)=w*(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^1;
        D1{face}(j,24)=w*(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^2;
        D1{face}(j,25)=w*(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^3;
        D1{face}(j,26)=w*(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^4;
        D1{face}(j,27)=w*(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^5;
        elseif extended_stencil && ny==0
        D1{face}(j,22)=w*(cells(cell,1)-faces(face,1))^6;
        D1{face}(j,23)=w*(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^1;
        D1{face}(j,24)=w*(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^2;
        D1{face}(j,25)=w*(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^3;
        D1{face}(j,26)=w*(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^4;
        D1{face}(j,27)=w*(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^5;
        end
        
        
        % Matriz D %
        D{face}(j,1)=1;
        %
        D{face}(j,2)=(cells(cell,1)-faces(face,1)+inteiro);
        D{face}(j,3)=(cells(cell,2)-faces(face,2)+inteiro);
        %
        D{face}(j,4)=(cells(cell,1)-faces(face,1)+inteiro)^2;
        D{face}(j,5)=(cells(cell,2)-faces(face,2)+inteiro)^2;
        D{face}(j,6)=(cells(cell,1)-faces(face,1)+inteiro)*(cells(cell,2)-faces(face,2)+inteiro);
        %
        D{face}(j,7)=(cells(cell,1)-faces(face,1)+inteiro)^3;
        D{face}(j,8)=(cells(cell,2)-faces(face,2)+inteiro)^3;
        D{face}(j,9)=(cells(cell,1)-faces(face,1)+inteiro)^2*(cells(cell,2)-faces(face,2)+inteiro);
        D{face}(j,10)=(cells(cell,1)-faces(face,1)+inteiro)^1*(cells(cell,2)-faces(face,2)+inteiro)^2;
        % 
        D{face}(j,11)=(cells(cell,1)-faces(face,1)+inteiro)^4;
        D{face}(j,12)=(cells(cell,2)-faces(face,2)+inteiro)^4;
        D{face}(j,13)=(cells(cell,1)-faces(face,1)+inteiro)^3*(cells(cell,2)-faces(face,2)+inteiro)^1;
        D{face}(j,14)=(cells(cell,1)-faces(face,1)+inteiro)^2*(cells(cell,2)-faces(face,2)+inteiro)^2;
        D{face}(j,15)=(cells(cell,1)-faces(face,1)+inteiro)^1*(cells(cell,2)-faces(face,2)+inteiro)^3;
        %
        D{face}(j,16)=(cells(cell,1)-faces(face,1)+inteiro)^5;
        D{face}(j,17)=(cells(cell,2)-faces(face,2)+inteiro)^5;
        D{face}(j,18)=(cells(cell,1)-faces(face,1)+inteiro)^4*(cells(cell,2)-faces(face,2)+inteiro)^1;
        D{face}(j,19)=(cells(cell,1)-faces(face,1)+inteiro)^3*(cells(cell,2)-faces(face,2)+inteiro)^2;
        D{face}(j,20)=(cells(cell,1)-faces(face,1)+inteiro)^2*(cells(cell,2)-faces(face,2)+inteiro)^3;
        D{face}(j,21)=(cells(cell,1)-faces(face,1)+inteiro)^1*(cells(cell,2)-faces(face,2)+inteiro)^4;
        
        
        if extended_stencil && nx==0       
        D{face}(j,22)=(cells(cell,2)-faces(face,2))^6;
        D{face}(j,23)=(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^1;
        D{face}(j,24)=(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^2;
        D{face}(j,25)=(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^3;
        D{face}(j,26)=(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^4;
        D{face}(j,27)=(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^5;
        elseif extended_stencil && ny==0
        D{face}(j,22)=(cells(cell,1)-faces(face,1))^6;
        D{face}(j,23)=(cells(cell,1)-faces(face,1))^5*(cells(cell,2)-faces(face,2))^1;
        D{face}(j,24)=(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^2;
        D{face}(j,25)=(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^3;
        D{face}(j,26)=(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^4;
        D{face}(j,27)=(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^5;
        end
        
        
        
    end
    %
    %% Construção da Matriz D para cada Face, Faces do Stencil %
    n_face=stencil_size(face,3);
    %
    if n_face~=0
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
%                     dx=0.0000001;
%                     dy=0.0000001;
                end
                w=1/(sqrt(dx^2+dy^2))^6;
            else
                w=1;
            end
            %
 %% NEUMAN     
 
 
            normal_x=face_bound(facestencil,2)*face_area(facestencil,1);
            normal_y=face_bound(facestencil,3)*face_area(facestencil,1);
            area_robin=1*face_area(facestencil,1);
%             if facestencil<cell_side+1
%             normal_x=-1*face_area(facestencil,1);    
%             end
            
            if ~dimensional_correction
                      normal_x=normal_x/face_area(facestencil,1);
                      normal_y=normal_y/face_area(facestencil,1);
                      area_robin=1;
              end 
            
            
            
            if neuman && face_bound(facestencil,4)==1
            D1{face}(n_cell+j,1)=0;
            %
            D1{face}(n_cell+j,2)=w*normal_x;
            D1{face}(n_cell+j,3)=w*normal_y;
            %
            D1{face}(n_cell+j,4)=2*w*(faces(facestencil,1)-faces(face,1))*normal_x;
            D1{face}(n_cell+j,5)=2*w*(faces(facestencil,2)-faces(face,2))*normal_y;
            D1{face}(n_cell+j,6)=w*(faces(facestencil,1)-faces(face,1))*normal_y+w*(faces(facestencil,2)-faces(face,2))*normal_x;
            %
            D1{face}(n_cell+j,7)=(3*w*(faces(facestencil,1)-faces(face,1))^2)*normal_x;
            D1{face}(n_cell+j,8)=(3*w*(faces(facestencil,2)-faces(face,2))^2)*normal_y;
            D1{face}(n_cell+j,9)=w*2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))*normal_x+w*((faces(facestencil,1)-faces(face,1))^2)*normal_y;
            D1{face}(n_cell+j,10)=w*((faces(facestencil,2)-faces(face,2))^2)*normal_x+w*2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))*normal_y;
            %
            D1{face}(n_cell+j,11)=4*w*(faces(facestencil,1)-faces(face,1))^3*normal_x;
            D1{face}(n_cell+j,12)=4*w*(faces(facestencil,2)-faces(face,2))^3*normal_y;           
            D1{face}(n_cell+j,13)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_x+w*(faces(facestencil,1)-faces(face,1))^3*normal_y;
            D1{face}(n_cell+j,14)=2*w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_y;
            D1{face}(n_cell+j,15)=w*(faces(facestencil,2)-faces(face,2))^3*normal_x+3*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_y;

            %
            D1{face}(n_cell+j,16)=5*w*(faces(facestencil,1)-faces(face,1))^4*normal_x;
            D1{face}(n_cell+j,17)=5*w*(faces(facestencil,2)-faces(face,2))^4*normal_y;
            D1{face}(n_cell+j,18)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1*normal_x+w*(faces(facestencil,1)-faces(face,1))^4*normal_y;
            D1{face}(n_cell+j,19)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1*normal_y;
            D1{face}(n_cell+j,20)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3*normal_x+3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2*normal_y;
            D1{face}(n_cell+j,21)=w*(faces(facestencil,2)-faces(face,2))^4*normal_x+4*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3*normal_y;
            %
            
            if extended_stencil && nx==0       
        D1{face}(n_cell+j,22)=0;
        D1{face}(n_cell+j,23)=5*w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,24)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,25)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D1{face}(n_cell+j,26)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D1{face}(n_cell+j,27)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        elseif extended_stencil && ny==0
        D1{face}(n_cell+j,22)=6*w*(faces(facestencil,1)-faces(face,1))^5*normal_x;
        D1{face}(n_cell+j,23)=5*w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,24)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,25)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D1{face}(n_cell+j,26)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D1{face}(n_cell+j,27)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        end
            
            
            
            % Matriz D %
           D{face}(n_cell+j,1)=0;
            %
            D{face}(n_cell+j,2)=normal_x;
            D{face}(n_cell+j,3)=normal_y;
            %
            D{face}(n_cell+j,4)=2*(faces(facestencil,1)-faces(face,1))*normal_x;
            D{face}(n_cell+j,5)=2*(faces(facestencil,2)-faces(face,2))*normal_y;
            D{face}(n_cell+j,6)=(faces(facestencil,1)-faces(face,1))*normal_y+(faces(facestencil,2)-faces(face,2))*normal_x;
            %
            D{face}(n_cell+j,7)=(3*(faces(facestencil,1)-faces(face,1))^2)*normal_x;
            D{face}(n_cell+j,8)=(3*(faces(facestencil,2)-faces(face,2))^2)*normal_y;
            D{face}(n_cell+j,9)=2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))*normal_x+((faces(facestencil,1)-faces(face,1))^2)*normal_y;
            D{face}(n_cell+j,10)=((faces(facestencil,2)-faces(face,2))^2)*normal_x+2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))*normal_y;
            %
            D{face}(n_cell+j,11)=4*(faces(facestencil,1)-faces(face,1))^3*normal_x;
            D{face}(n_cell+j,12)=4*(faces(facestencil,2)-faces(face,2))^3*normal_y;           
            D{face}(n_cell+j,13)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_x+(faces(facestencil,1)-faces(face,1))^3*normal_y;
            D{face}(n_cell+j,14)=2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_y;
            D{face}(n_cell+j,15)=(faces(facestencil,2)-faces(face,2))^3*normal_x+3*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_y;

            %
            D{face}(n_cell+j,16)=5*(faces(facestencil,1)-faces(face,1))^4*normal_x;
            D{face}(n_cell+j,17)=5*(faces(facestencil,2)-faces(face,2))^4*normal_y;
            D{face}(n_cell+j,18)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1*normal_x+(faces(facestencil,1)-faces(face,1))^4*normal_y;
            D{face}(n_cell+j,19)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1*normal_y;
            D{face}(n_cell+j,20)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3*normal_x+3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2*normal_y;
            D{face}(n_cell+j,21)=(faces(facestencil,2)-faces(face,2))^4*normal_x+4*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3*normal_y;
            
            
        if extended_stencil && nx==0       
        D{face}(n_cell+j,22)=0;
        D{face}(n_cell+j,23)=5*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,24)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,25)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D{face}(n_cell+j,26)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D{face}(n_cell+j,27)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        elseif extended_stencil && ny==0
        D{face}(n_cell+j,22)=6*(faces(facestencil,1)-faces(face,1))^5*normal_x;
        D{face}(n_cell+j,23)=5*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,24)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,25)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D{face}(n_cell+j,26)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D{face}(n_cell+j,27)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        end    
                
%% ROBIN


            elseif robin && face_bound(facestencil,4)==2
            D1{face}(n_cell+j,1)=0;
            %
            D1{face}(n_cell+j,2)=w*normal_x;
            D1{face}(n_cell+j,3)=w*normal_y;
            %
            D1{face}(n_cell+j,4)=2*w*(faces(facestencil,1)-faces(face,1))*normal_x;
            D1{face}(n_cell+j,5)=2*w*(faces(facestencil,2)-faces(face,2))*normal_y;
            D1{face}(n_cell+j,6)=w*(faces(facestencil,1)-faces(face,1))*normal_y+w*(faces(facestencil,2)-faces(face,2))*normal_x;
            %
            D1{face}(n_cell+j,7)=(3*w*(faces(facestencil,1)-faces(face,1))^2)*normal_x;
            D1{face}(n_cell+j,8)=(3*w*(faces(facestencil,2)-faces(face,2))^2)*normal_y;
            D1{face}(n_cell+j,9)=w*2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))*normal_x+w*((faces(facestencil,1)-faces(face,1))^2)*normal_y;
            D1{face}(n_cell+j,10)=w*((faces(facestencil,2)-faces(face,2))^2)*normal_x+w*2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))*normal_y;
            %
            D1{face}(n_cell+j,11)=4*w*(faces(facestencil,1)-faces(face,1))^3*normal_x;
            D1{face}(n_cell+j,12)=4*w*(faces(facestencil,2)-faces(face,2))^3*normal_y;           
            D1{face}(n_cell+j,13)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_x+w*(faces(facestencil,1)-faces(face,1))^3*normal_y;
            D1{face}(n_cell+j,14)=2*w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_y;
            D1{face}(n_cell+j,15)=w*(faces(facestencil,2)-faces(face,2))^3*normal_x+3*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_y;

            %
            D1{face}(n_cell+j,16)=5*w*(faces(facestencil,1)-faces(face,1))^4*normal_x;
            D1{face}(n_cell+j,17)=5*w*(faces(facestencil,2)-faces(face,2))^4*normal_y;
            D1{face}(n_cell+j,18)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1*normal_x+w*(faces(facestencil,1)-faces(face,1))^4*normal_y;
            D1{face}(n_cell+j,19)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1*normal_y;
            D1{face}(n_cell+j,20)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3*normal_x+3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2*normal_y;
            D1{face}(n_cell+j,21)=w*(faces(facestencil,2)-faces(face,2))^4*normal_x+4*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3*normal_y;
            %
            
            if extended_stencil && nx==0       
        D1{face}(n_cell+j,22)=0;
        D1{face}(n_cell+j,23)=5*w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,24)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,25)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D1{face}(n_cell+j,26)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D1{face}(n_cell+j,27)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        elseif extended_stencil && ny==0
        D1{face}(n_cell+j,22)=6*w*(faces(facestencil,1)-faces(face,1))^5*normal_x;
        D1{face}(n_cell+j,23)=5*w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D1{face}(n_cell+j,24)=4*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D1{face}(n_cell+j,25)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D1{face}(n_cell+j,26)=2*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D1{face}(n_cell+j,27)=1*w*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        end
            
            
            
            % Matriz D %
           D{face}(n_cell+j,1)=0;
            %
            D{face}(n_cell+j,2)=normal_x;
            D{face}(n_cell+j,3)=normal_y;
            %
            D{face}(n_cell+j,4)=2*(faces(facestencil,1)-faces(face,1))*normal_x;
            D{face}(n_cell+j,5)=2*(faces(facestencil,2)-faces(face,2))*normal_y;
            D{face}(n_cell+j,6)=(faces(facestencil,1)-faces(face,1))*normal_y+(faces(facestencil,2)-faces(face,2))*normal_x;
            %
            D{face}(n_cell+j,7)=(3*(faces(facestencil,1)-faces(face,1))^2)*normal_x;
            D{face}(n_cell+j,8)=(3*(faces(facestencil,2)-faces(face,2))^2)*normal_y;
            D{face}(n_cell+j,9)=2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))*normal_x+((faces(facestencil,1)-faces(face,1))^2)*normal_y;
            D{face}(n_cell+j,10)=((faces(facestencil,2)-faces(face,2))^2)*normal_x+2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))*normal_y;
            %
            D{face}(n_cell+j,11)=4*(faces(facestencil,1)-faces(face,1))^3*normal_x;
            D{face}(n_cell+j,12)=4*(faces(facestencil,2)-faces(face,2))^3*normal_y;           
            D{face}(n_cell+j,13)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_x+(faces(facestencil,1)-faces(face,1))^3*normal_y;
            D{face}(n_cell+j,14)=2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_y;
            D{face}(n_cell+j,15)=(faces(facestencil,2)-faces(face,2))^3*normal_x+3*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_y;

            %
            D{face}(n_cell+j,16)=5*(faces(facestencil,1)-faces(face,1))^4*normal_x;
            D{face}(n_cell+j,17)=5*(faces(facestencil,2)-faces(face,2))^4*normal_y;
            D{face}(n_cell+j,18)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1*normal_x+(faces(facestencil,1)-faces(face,1))^4*normal_y;
            D{face}(n_cell+j,19)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1*normal_y;
            D{face}(n_cell+j,20)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3*normal_x+3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2*normal_y;
            D{face}(n_cell+j,21)=(faces(facestencil,2)-faces(face,2))^4*normal_x+4*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3*normal_y;
            
            
        if extended_stencil && nx==0       
        D{face}(n_cell+j,22)=0;
        D{face}(n_cell+j,23)=5*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,24)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,25)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D{face}(n_cell+j,26)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D{face}(n_cell+j,27)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        elseif extended_stencil && ny==0
        D{face}(n_cell+j,22)=6*(faces(facestencil,1)-faces(face,1))^5*normal_x;
        D{face}(n_cell+j,23)=5*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1*normal_x;
        D{face}(n_cell+j,24)=4*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2*normal_x;
        D{face}(n_cell+j,25)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3*normal_x;
        D{face}(n_cell+j,26)=2*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4*normal_x;
        D{face}(n_cell+j,27)=1*(faces(facestencil,1)-faces(face,1))^0*(faces(facestencil,2)-faces(face,2))^5*normal_x;
        end 
%%%%%%%%%%
            D1{face}(n_cell+j,1)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,1)+ area_robin*w*1;
            %
            D1{face}(n_cell+j,2)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,2)+area_robin*w*(faces(facestencil,1)-faces(face,1));
            D1{face}(n_cell+j,3)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,3)+area_robin*w*(faces(facestencil,2)-faces(face,2));
            %
            D1{face}(n_cell+j,4)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,4)+area_robin*w*(faces(facestencil,1)-faces(face,1))^2;
            D1{face}(n_cell+j,5)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,5)+area_robin*w*(faces(facestencil,2)-faces(face,2))^2;
            D1{face}(n_cell+j,6)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,6)+area_robin*w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2));
            %
            D1{face}(n_cell+j,7)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,7)+area_robin*w*(faces(facestencil,1)-faces(face,1))^3;
            D1{face}(n_cell+j,8)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,8)+area_robin*w*(faces(facestencil,2)-faces(face,2))^3;
            D1{face}(n_cell+j,9)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,9)+area_robin*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2));
            D1{face}(n_cell+j,10)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,10)+area_robin*w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2;
            %
            D1{face}(n_cell+j,11)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,11)+area_robin*w*(faces(facestencil,1)-faces(face,1))^4;
            D1{face}(n_cell+j,12)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,12)+area_robin*w*(faces(facestencil,2)-faces(face,2))^4;
            D1{face}(n_cell+j,13)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,13)+area_robin*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1;
            D1{face}(n_cell+j,14)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,14)+area_robin*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2;
            D1{face}(n_cell+j,15)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,15)+area_robin*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3;
            %
            D1{face}(n_cell+j,16)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,16)+area_robin*w*(faces(facestencil,1)-faces(face,1))^5;
            D1{face}(n_cell+j,17)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,17)+area_robin*w*(faces(facestencil,2)-faces(face,2))^5;
            D1{face}(n_cell+j,18)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,18)+area_robin*w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1;
            D1{face}(n_cell+j,19)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,19)+area_robin*w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2;
            D1{face}(n_cell+j,20)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,20)+area_robin*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3;
            D1{face}(n_cell+j,21)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,21)+area_robin*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4;
            %
            
            if extended_stencil && nx==0       
        D1{face}(n_cell+j,22)=w*(faces(facestencil,2)-faces(face,2))^6;
        D1{face}(n_cell+j,23)=w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,24)=w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,25)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,26)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4;
        D1{face}(n_cell+j,27)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5;
        elseif extended_stencil && ny==0
        D1{face}(n_cell+j,22)=w*(faces(facestencil,1)-faces(face,1))^6;
        D1{face}(n_cell+j,23)=w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,24)=w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,25)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,26)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4;
        D1{face}(n_cell+j,27)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5;
        end
            
            
            
            % Matriz D %
            D{face}(n_cell+j,1)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,1)+area_robin*1;
            %
            D{face}(n_cell+j,2)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,2)+area_robin*(faces(facestencil,1)-faces(face,1));
            D{face}(n_cell+j,3)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,3)+area_robin*(faces(facestencil,2)-faces(face,2));
            %
            D{face}(n_cell+j,4)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,4)+area_robin*(faces(facestencil,1)-faces(face,1))^2;
            D{face}(n_cell+j,5)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,5)+area_robin*(faces(facestencil,2)-faces(face,2))^2;
            D{face}(n_cell+j,6)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,6)+area_robin*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2));
            %
            D{face}(n_cell+j,7)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,7)+area_robin*(faces(facestencil,1)-faces(face,1))^3;
            D{face}(n_cell+j,8)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,8)+area_robin*(faces(facestencil,2)-faces(face,2))^3;
            D{face}(n_cell+j,9)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,9)+area_robin*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2));
            D{face}(n_cell+j,10)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,10)+area_robin*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2;
            %
            D{face}(n_cell+j,11)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,11)+area_robin*(faces(facestencil,1)-faces(face,1))^4;
            D{face}(n_cell+j,12)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,12)+area_robin*(faces(facestencil,2)-faces(face,2))^4;
            D{face}(n_cell+j,13)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,13)+area_robin*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1;
            D{face}(n_cell+j,14)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,14)+area_robin*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2;
            D{face}(n_cell+j,15)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,15)+area_robin*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3;
            %
            D{face}(n_cell+j,16)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,16)+area_robin*(faces(facestencil,1)-faces(face,1))^5;
            D{face}(n_cell+j,17)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,17)+area_robin*(faces(facestencil,2)-faces(face,2))^5;
            D{face}(n_cell+j,18)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,18)+area_robin*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^1;
            D{face}(n_cell+j,19)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,19)+area_robin*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^2;
            D{face}(n_cell+j,20)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,20)+area_robin*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^3;
            D{face}(n_cell+j,21)= (gamma_diff/u_convec_x)*D{face}(n_cell+j,21)+area_robin*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^4;
            
            
        if extended_stencil && nx==0       
        D{face}(n_cell+j,22)=(faces(facestencil,2)-faces(face,2))^6;
        D{face}(n_cell+j,23)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,24)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,25)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,26)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,27)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5;
        elseif extended_stencil && ny==0
        D{face}(n_cell+j,22)=(faces(facestencil,1)-faces(face,1))^6;
        D{face}(n_cell+j,23)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,24)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,25)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,26)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,27)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5;
        end





%% Drichlet
            else
            % Matriz D tendo em conta a ponderação %
            D1{face}(n_cell+j,1)=w*1;
            %
            D1{face}(n_cell+j,2)=w*(faces(facestencil,1)-faces(face,1)+inteiro);
            D1{face}(n_cell+j,3)=w*(faces(facestencil,2)-faces(face,2)+inteiro);
            %
            D1{face}(n_cell+j,4)=w*(faces(facestencil,1)-faces(face,1)+inteiro)^2;
            D1{face}(n_cell+j,5)=w*(faces(facestencil,2)-faces(face,2)+inteiro)^2;
            D1{face}(n_cell+j,6)=w*(faces(facestencil,1)-faces(face,1)+inteiro)*(faces(facestencil,2)-faces(face,2)+inteiro);
            %
            D1{face}(n_cell+j,7)=w*(faces(facestencil,1)-faces(face,1)+inteiro)^3;
            D1{face}(n_cell+j,8)=w*(faces(facestencil,2)-faces(face,2)+inteiro)^3;
            D1{face}(n_cell+j,9)=w*(faces(facestencil,1)-faces(face,1)+inteiro)^2*(faces(facestencil,2)-faces(face,2)+inteiro);
            D1{face}(n_cell+j,10)=w*(faces(facestencil,1)-faces(face,1)+inteiro)*(faces(facestencil,2)-faces(face,2)+inteiro)^2;
            %
            D1{face}(n_cell+j,11)=w*(faces(facestencil,1)-faces(face,1)+inteiro)^4;
            D1{face}(n_cell+j,12)=w*(faces(facestencil,2)-faces(face,2)+inteiro)^4;
            D1{face}(n_cell+j,13)=w*(faces(facestencil,1)-faces(face,1)+inteiro)^3*(faces(facestencil,2)-faces(face,2)+inteiro)^1;
            D1{face}(n_cell+j,14)=w*(faces(facestencil,1)-faces(face,1)+inteiro)^2*(faces(facestencil,2)-faces(face,2)+inteiro)^2;
            D1{face}(n_cell+j,15)=w*(faces(facestencil,1)-faces(face,1)+inteiro)^1*(faces(facestencil,2)-faces(face,2)+inteiro)^3;
            %
            D1{face}(n_cell+j,16)=w*(faces(facestencil,1)-faces(face,1)+inteiro)^5;
            D1{face}(n_cell+j,17)=w*(faces(facestencil,2)-faces(face,2)+inteiro)^5;
            D1{face}(n_cell+j,18)=w*(faces(facestencil,1)-faces(face,1)+inteiro)^4*(faces(facestencil,2)-faces(face,2)+inteiro)^1;
            D1{face}(n_cell+j,19)=w*(faces(facestencil,1)-faces(face,1)+inteiro)^3*(faces(facestencil,2)-faces(face,2)+inteiro)^2;
            D1{face}(n_cell+j,20)=w*(faces(facestencil,1)-faces(face,1)+inteiro)^2*(faces(facestencil,2)-faces(face,2)+inteiro)^3;
            D1{face}(n_cell+j,21)=w*(faces(facestencil,1)-faces(face,1)+inteiro)^1*(faces(facestencil,2)-faces(face,2)+inteiro)^4;
            %
            
            if extended_stencil && nx==0       
        D1{face}(n_cell+j,22)=w*(faces(facestencil,2)-faces(face,2))^6;
        D1{face}(n_cell+j,23)=w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,24)=w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,25)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,26)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4;
        D1{face}(n_cell+j,27)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5;
        elseif extended_stencil && ny==0
        D1{face}(n_cell+j,22)=w*(faces(facestencil,1)-faces(face,1))^6;
        D1{face}(n_cell+j,23)=w*(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1;
        D1{face}(n_cell+j,24)=w*(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2;
        D1{face}(n_cell+j,25)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3;
        D1{face}(n_cell+j,26)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4;
        D1{face}(n_cell+j,27)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5;
        end
            
            
            
            % Matriz D %
            D{face}(n_cell+j,1)=1;
            %
            D{face}(n_cell+j,2)=(faces(facestencil,1)-faces(face,1)+inteiro);
            D{face}(n_cell+j,3)=(faces(facestencil,2)-faces(face,2)+inteiro);
            %
            D{face}(n_cell+j,4)=(faces(facestencil,1)-faces(face,1)+inteiro)^2;
            D{face}(n_cell+j,5)=(faces(facestencil,2)-faces(face,2)+inteiro)^2;
            D{face}(n_cell+j,6)=(faces(facestencil,1)-faces(face,1)+inteiro)*(faces(facestencil,2)-faces(face,2)+inteiro);
            %
            D{face}(n_cell+j,7)=(faces(facestencil,1)-faces(face,1)+inteiro)^3;
            D{face}(n_cell+j,8)=(faces(facestencil,2)-faces(face,2)+inteiro)^3;
            D{face}(n_cell+j,9)=(faces(facestencil,1)-faces(face,1)+inteiro)^2*(faces(facestencil,2)-faces(face,2)+inteiro);
            D{face}(n_cell+j,10)=(faces(facestencil,1)-faces(face,1)+inteiro)*(faces(facestencil,2)-faces(face,2)+inteiro)^2;
            %
            D{face}(n_cell+j,11)=(faces(facestencil,1)-faces(face,1)+inteiro)^4;
            D{face}(n_cell+j,12)=(faces(facestencil,2)-faces(face,2)+inteiro)^4;
            D{face}(n_cell+j,13)=(faces(facestencil,1)-faces(face,1)+inteiro)^3*(faces(facestencil,2)-faces(face,2)+inteiro)^1;
            D{face}(n_cell+j,14)=(faces(facestencil,1)-faces(face,1)+inteiro)^2*(faces(facestencil,2)-faces(face,2)+inteiro)^2;
            D{face}(n_cell+j,15)=(faces(facestencil,1)-faces(face,1)+inteiro)^1*(faces(facestencil,2)-faces(face,2)+inteiro)^3;
            %
            D{face}(n_cell+j,16)=(faces(facestencil,1)-faces(face,1)+inteiro)^5;
            D{face}(n_cell+j,17)=(faces(facestencil,2)-faces(face,2)+inteiro)^5;
            D{face}(n_cell+j,18)=(faces(facestencil,1)-faces(face,1)+inteiro)^4*(faces(facestencil,2)-faces(face,2)+inteiro)^1;
            D{face}(n_cell+j,19)=(faces(facestencil,1)-faces(face,1)+inteiro)^3*(faces(facestencil,2)-faces(face,2)+inteiro)^2;
            D{face}(n_cell+j,20)=(faces(facestencil,1)-faces(face,1)+inteiro)^2*(faces(facestencil,2)-faces(face,2)+inteiro)^3;
            D{face}(n_cell+j,21)=(faces(facestencil,1)-faces(face,1)+inteiro)^1*(faces(facestencil,2)-faces(face,2)+inteiro)^4;
            
            
        if extended_stencil && nx==0       
        D{face}(n_cell+j,22)=(faces(facestencil,2)-faces(face,2))^6;
        D{face}(n_cell+j,23)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,24)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,25)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,26)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,27)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5;
        elseif extended_stencil && ny==0
        D{face}(n_cell+j,22)=(faces(facestencil,1)-faces(face,1))^6;
        D{face}(n_cell+j,23)=(faces(facestencil,1)-faces(face,1))^5*(faces(facestencil,2)-faces(face,2))^1;
        D{face}(n_cell+j,24)=(faces(facestencil,1)-faces(face,1))^4*(faces(facestencil,2)-faces(face,2))^2;
        D{face}(n_cell+j,25)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^3;
        D{face}(n_cell+j,26)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^4;
        D{face}(n_cell+j,27)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^5;
        end
            
            end 
        end
    end
    %
    %% Determinação da Matriz T %%
    % D*c=phi -> c=inv(D'D)*D'*phi -> T=inv(D'D)*D' %
    T{face}=inv(D{face}'*D1{face})*D1{face}';
end
pause(0.00001)
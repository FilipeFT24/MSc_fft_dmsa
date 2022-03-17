function [T,D]=Reconstruction4thOrder(stencil_cells,stencil_faces,stencil_size,ponderado)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                  06 de Outubro de 2016                                            %
%                                                                                                   %
% Função que implementa a Solução Numérica a partir do Minimos Quadrados de 2ª Ordem                %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global L Lref cell_side vert_side cell_num face_num vert_num w;
global verts cells cell_verts cell_faces cell_vol faces face_area face_vert cell_norm;
global face_bound cell_bound face_cells vert_cells;
global phi lap_phi  phi_faces flux_phi_faces face_w_gauss extended_stencil neuman dimensional_correction robin u_convec_x u_convec_y gamma_diff;
%
for i=1:face_num
    %
    face=i;
    %
    % Normal Exterior da Face %
    vert=face_vert(face,:);
    nx=verts(vert(1),1)-verts(vert(2),1);
    ny=verts(vert(1),2)-verts(vert(2),2);
    if nx~=0
        nx=1;
    elseif ny~=0
        ny=1;
    end

            
    %% Construção da Matriz D para cada Face, Celulas do Stencil %%
    n_cell=stencil_size(face,2);
    %
    for j=1:n_cell
        celula=stencil_cells(face,j);
        %
        if ponderado==true
            dx=(cells(celula,1)-faces(face,1));
            dy=(cells(celula,2)-faces(face,2));
            w=1/(sqrt(dx^2+dy^2))^4;
        else
            w=1;
        end
        %
        % Matriz D tendo em conta a ponderação %
        D1{face}(j,1)=w*1;
        %
        D1{face}(j,2)=w*(cells(celula,1)-faces(face,1));
        D1{face}(j,3)=w*(cells(celula,2)-faces(face,2));
        %
        D1{face}(j,4)=w*(cells(celula,1)-faces(face,1))^2;
        D1{face}(j,5)=w*(cells(celula,2)-faces(face,2))^2;
        D1{face}(j,6)=w*(cells(celula,1)-faces(face,1))*(cells(celula,2)-faces(face,2));
        %
        D1{face}(j,7)=w*(cells(celula,1)-faces(face,1))^3;
        D1{face}(j,8)=w*(cells(celula,2)-faces(face,2))^3;
        D1{face}(j,9)=w*(cells(celula,1)-faces(face,1))^2*(cells(celula,2)-faces(face,2));
        D1{face}(j,10)=w*(cells(celula,1)-faces(face,1))*(cells(celula,2)-faces(face,2))^2;
        
        hold on;
        plot(cells(celula,1),cells(celula,2),'ob');
        
        
        if extended_stencil && nx==0
        D1{face}(j,11)=w*(cells(celula,2)-faces(face,2))^4;
        D1{face}(j,12)=w*(cells(celula,1)-faces(face,1))^3*(cells(celula,2)-faces(face,2));
        D1{face}(j,13)=w*(cells(celula,1)-faces(face,1))^2*(cells(celula,2)-faces(face,2))^2;
        D1{face}(j,14)=w*(cells(celula,1)-faces(face,1))^1*(cells(celula,2)-faces(face,2))^3;
        elseif extended_stencil && ny==0
        D1{face}(j,11)=w*(cells(celula,1)-faces(face,1))^4;
        D1{face}(j,12)=w*(cells(celula,1)-faces(face,1))^3*(cells(celula,2)-faces(face,2));
        D1{face}(j,13)=w*(cells(celula,1)-faces(face,1))^2*(cells(celula,2)-faces(face,2))^2;
        D1{face}(j,14)=w*(cells(celula,1)-faces(face,1))^1*(cells(celula,2)-faces(face,2))^3;
        end
        
        % Matriz D %
        D{face}(j,1)=1;
        %
        D{face}(j,2)=(cells(celula,1)-faces(face,1));
        D{face}(j,3)=(cells(celula,2)-faces(face,2));
        %
        D{face}(j,4)=(cells(celula,1)-faces(face,1))^2;
        D{face}(j,5)=(cells(celula,2)-faces(face,2))^2;
        D{face}(j,6)=(cells(celula,1)-faces(face,1))*(cells(celula,2)-faces(face,2));
        %
        D{face}(j,7)=(cells(celula,1)-faces(face,1))^3;
        D{face}(j,8)=(cells(celula,2)-faces(face,2))^3;
        D{face}(j,9)=(cells(celula,1)-faces(face,1))^2*(cells(celula,2)-faces(face,2));
        D{face}(j,10)=(cells(celula,1)-faces(face,1))*(cells(celula,2)-faces(face,2))^2;
        
        
       if extended_stencil && nx==0
            
    
        D{face}(j,11)=(cells(celula,2)-faces(face,2))^4;
        D{face}(j,12)=(cells(celula,1)-faces(face,1))^3*(cells(celula,2)-faces(face,2));
        D{face}(j,13)=(cells(celula,1)-faces(face,1))^2*(cells(celula,2)-faces(face,2))^2;
        D{face}(j,14)=(cells(celula,1)-faces(face,1))^1*(cells(celula,2)-faces(face,2))^3;
        
        elseif extended_stencil && ny==0
           
        D{face}(j,11)=(cells(celula,1)-faces(face,1))^4;
        D{face}(j,12)=(cells(celula,1)-faces(face,1))^3*(cells(celula,2)-faces(face,2));
        D{face}(j,13)=(cells(celula,1)-faces(face,1))^2*(cells(celula,2)-faces(face,2))^2;
        D{face}(j,14)=(cells(celula,1)-faces(face,1))^1*(cells(celula,2)-faces(face,2))^3;
        
        end
        
        
    end
    %
    %% Construção da Matriz D para cada Face, Faces do Stencil %%
    n_face=stencil_size(face,3);
    %
    if n_face~=0
        for j=1:n_face
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
                w=1/(sqrt(dx^2+dy^2))^4;
            else
                w=1;
            end
            
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
  if normal_y~=0 
      pause()
  end
        
                
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
        D1{face}(n_cell+j,9)=2*w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))*normal_x;
        D1{face}(n_cell+j,10)=w*(faces(facestencil,2)-faces(face,2))^2*normal_x;
            %
            
        if extended_stencil && nx==0
            D1{face}(n_cell+j,11)=4*w*(faces(facestencil,2)-faces(face,2))^3*normal_y;
            D1{face}(n_cell+j,12)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_x+w*(faces(facestencil,1)-faces(face,1))^3*normal_y;
            D1{face}(n_cell+j,13)=2*w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_y;
            D1{face}(n_cell+j,14)=w*(faces(facestencil,2)-faces(face,2))^3*normal_x+3*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_y;
        elseif extended_stencil && ny==0
            D1{face}(n_cell+j,11)=4*w*(faces(facestencil,1)-faces(face,1))^3*normal_x;
            D1{face}(n_cell+j,12)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_x+w*(faces(facestencil,1)-faces(face,1))^3*normal_y;
            D1{face}(n_cell+j,13)=2*w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_y;
            D1{face}(n_cell+j,14)=w*(faces(facestencil,2)-faces(face,2))^3*normal_x+3*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_y;
        end
            
            
             % Matriz D %
        D{face}(n_cell+j,1)=0;
        %
        D{face}(n_cell+j,2)=1*normal_x;
        D{face}(n_cell+j,3)=0;
        %
        D{face}(n_cell+j,4)=2*(faces(facestencil,1)-faces(face,1))*normal_x;
        D{face}(n_cell+j,5)=0;
        D{face}(n_cell+j,6)=(faces(facestencil,2)-faces(face,2))*normal_x;
        %
        D{face}(n_cell+j,7)=3*(faces(facestencil,1)-faces(face,1))^2*normal_x;
        D{face}(n_cell+j,8)=0;
        D{face}(n_cell+j,9)=2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))*normal_x;
        D{face}(n_cell+j,10)=(faces(facestencil,2)-faces(face,2))^2*normal_x;
            
           if extended_stencil && nx==0
            D{face}(n_cell+j,11)=4*(faces(facestencil,2)-faces(face,2))^3*normal_y;
            D{face}(n_cell+j,12)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_x+(faces(facestencil,1)-faces(face,1))^3*normal_y;
            D{face}(n_cell+j,13)=2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_y;
            D{face}(n_cell+j,14)=(faces(facestencil,2)-faces(face,2))^3*normal_x+3*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_y;
        elseif extended_stencil && ny==0
            D{face}(n_cell+j,11)=4*(faces(facestencil,1)-faces(face,1))^3*normal_x;
            D{face}(n_cell+j,12)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_x+(faces(facestencil,1)-faces(face,1))^3*normal_y;
            D{face}(n_cell+j,13)=2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_y;
            D{face}(n_cell+j,14)=(faces(facestencil,2)-faces(face,2))^3*normal_x+3*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_y;
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
        D1{face}(n_cell+j,9)=2*w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))*normal_x;
        D1{face}(n_cell+j,10)=w*(faces(facestencil,2)-faces(face,2))^2*normal_x;
            %
            
        if extended_stencil && nx==0
            D1{face}(n_cell+j,11)=4*w*(faces(facestencil,2)-faces(face,2))^3*normal_y;
            D1{face}(n_cell+j,12)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_x+w*(faces(facestencil,1)-faces(face,1))^3*normal_y;
            D1{face}(n_cell+j,13)=2*w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_y;
            D1{face}(n_cell+j,14)=w*(faces(facestencil,2)-faces(face,2))^3*normal_x+3*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_y;
        elseif extended_stencil && ny==0
            D1{face}(n_cell+j,11)=4*w*(faces(facestencil,1)-faces(face,1))^3*normal_x;
            D1{face}(n_cell+j,12)=3*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_x+w*(faces(facestencil,1)-faces(face,1))^3*normal_y;
            D1{face}(n_cell+j,13)=2*w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_y;
            D1{face}(n_cell+j,14)=w*(faces(facestencil,2)-faces(face,2))^3*normal_x+3*w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_y;
        end
            
            
             % Matriz D %
        D{face}(n_cell+j,1)=0;
        %
        D{face}(n_cell+j,2)=1*normal_x;
        D{face}(n_cell+j,3)=0;
        %
        D{face}(n_cell+j,4)=2*(faces(facestencil,1)-faces(face,1))*normal_x;
        D{face}(n_cell+j,5)=0;
        D{face}(n_cell+j,6)=(faces(facestencil,2)-faces(face,2))*normal_x;
        %
        D{face}(n_cell+j,7)=3*(faces(facestencil,1)-faces(face,1))^2*normal_x;
        D{face}(n_cell+j,8)=0;
        D{face}(n_cell+j,9)=2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))*normal_x;
        D{face}(n_cell+j,10)=(faces(facestencil,2)-faces(face,2))^2*normal_x;
            
           if extended_stencil && nx==0
            D{face}(n_cell+j,11)=4*(faces(facestencil,2)-faces(face,2))^3*normal_y;
            D{face}(n_cell+j,12)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_x+(faces(facestencil,1)-faces(face,1))^3*normal_y;
            D{face}(n_cell+j,13)=2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_y;
            D{face}(n_cell+j,14)=(faces(facestencil,2)-faces(face,2))^3*normal_x+3*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_y;
        elseif extended_stencil && ny==0
            D{face}(n_cell+j,11)=4*(faces(facestencil,1)-faces(face,1))^3*normal_x;
            D{face}(n_cell+j,12)=3*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_x+(faces(facestencil,1)-faces(face,1))^3*normal_y;
            D{face}(n_cell+j,13)=2*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2))^2*normal_x+2*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))*normal_y;
            D{face}(n_cell+j,14)=(faces(facestencil,2)-faces(face,2))^3*normal_x+3*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^2*normal_y;
           end      
 
        
           
 
           
             D1{face}(n_cell+j,1)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,1)+area_robin*w*1;
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
            
        if extended_stencil && nx==0
            D1{face}(n_cell+j,11)=w*(faces(facestencil,2)-faces(face,2))^4;
            D1{face}(n_cell+j,12)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1;
            D1{face}(n_cell+j,13)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2;
            D1{face}(n_cell+j,14)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3;
        elseif extended_stencil && ny==0
            D1{face}(n_cell+j,11)=w*(faces(facestencil,1)-faces(face,1))^4;
            D1{face}(n_cell+j,12)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1;
            D1{face}(n_cell+j,13)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2;
            D1{face}(n_cell+j,14)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3;
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
            
            
         if extended_stencil && nx==0
            D{face}(n_cell+j,11)=(faces(facestencil,2)-faces(face,2))^4;
            D{face}(n_cell+j,12)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1;
            D{face}(n_cell+j,13)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2;
            D{face}(n_cell+j,14)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3;
        elseif extended_stencil && ny==0
            D{face}(n_cell+j,11)=(faces(facestencil,1)-faces(face,1))^4;
            D{face}(n_cell+j,12)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1;
            D{face}(n_cell+j,13)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2;
            D{face}(n_cell+j,14)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3;
        end
            
                
 %% DRICHLET               
            else
            %
            % Matriz D tendo em conta a ponderação %
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
            
        if extended_stencil && nx==0
            D1{face}(n_cell+j,11)=w*(faces(facestencil,2)-faces(face,2))^4;
            D1{face}(n_cell+j,12)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1;
            D1{face}(n_cell+j,13)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2;
            D1{face}(n_cell+j,14)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3;
        elseif extended_stencil && ny==0
            D1{face}(n_cell+j,11)=w*(faces(facestencil,1)-faces(face,1))^4;
            D1{face}(n_cell+j,12)=w*(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1;
            D1{face}(n_cell+j,13)=w*(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2;
            D1{face}(n_cell+j,14)=w*(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3;
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
            
            
         if extended_stencil && nx==0
            D{face}(n_cell+j,11)=(faces(facestencil,2)-faces(face,2))^4;
            D{face}(n_cell+j,12)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1;
            D{face}(n_cell+j,13)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2;
            D{face}(n_cell+j,14)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3;
        elseif extended_stencil && ny==0
            D{face}(n_cell+j,11)=(faces(facestencil,1)-faces(face,1))^4;
            D{face}(n_cell+j,12)=(faces(facestencil,1)-faces(face,1))^3*(faces(facestencil,2)-faces(face,2))^1;
            D{face}(n_cell+j,13)=(faces(facestencil,1)-faces(face,1))^2*(faces(facestencil,2)-faces(face,2))^2;
            D{face}(n_cell+j,14)=(faces(facestencil,1)-faces(face,1))^1*(faces(facestencil,2)-faces(face,2))^3;
        end
            
            end
        end
    end
    %
    %% Determinação da Matriz T %%
    %
%     D*c=phi -> c=inv(D'D)*D'*phi -> T=inv(D'D)*D' %

    
T{face}=inv(D{face}'*D1{face})*D1{face}';
condf{face} = cond(D{face});

end

 pause(0.01)
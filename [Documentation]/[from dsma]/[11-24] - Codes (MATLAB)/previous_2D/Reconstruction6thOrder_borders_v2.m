function [T_border,D, constrained_source]=Reconstruction6thOrder_borders_v2(stencil_cells,stencil_faces,stencil_size,ponderado)
warning('on', 'all')
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
global face_bound cell_bound face_cells vert_cells dimensional_correction robin u_convec_x u_convec_y gamma_diff;
global phi lap_phi  phi_faces flux_phi_faces face_w_gauss G restos solution extended_stencil neuman;
%
contador535151=0;
q=1;

if face_w_gauss
for i=1:face_num
    %
    face=i;
          % Normal Exterior da Face %
    vert=face_vert(face,:);
    nx=verts(vert(1),1)-verts(vert(2),1);
    ny=verts(vert(1),2)-verts(vert(2),2);
    
    if face_bound(face,1)==1
    %% Pontos Gauss
    x(1)=G{face}(1,1);
    y(1)=G{face}(1,2);
                    %
                % Ponto de Gauss 2 %
                x(2)=G{face}(2,1);
                y(2)=G{face}(2,2);
 
                %
                % Ponto de Gauss 3 %
                x(3)=G{face}(3,1);
                y(3)=G{face}(3,2);
           
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
        soluti{face}(j) = SolutionDiffusion(solution,cells(cell,1),cells(cell,2),'anal');       
        % Matriz D tendo em conta a ponderação %
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
        D{face}(j,10)=(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^2;
        % 
        D{face}(j,11)=(cells(cell,1)-faces(face,1))^4;
        D{face}(j,12)=(cells(cell,2)-faces(face,2))^4;
        D{face}(j,13)=(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^1;
        D{face}(j,14)=(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^2;
        D{face}(j,15)=(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^3;
        %
        D{face}(j,16)=(cells(cell,1)-faces(face,1))^5;
        D{face}(j,17)=(cells(cell,2)-faces(face,2))^5;
        D{face}(j,18)=(cells(cell,1)-faces(face,1))^4*(cells(cell,2)-faces(face,2))^1;
        D{face}(j,19)=(cells(cell,1)-faces(face,1))^3*(cells(cell,2)-faces(face,2))^2;
        D{face}(j,20)=(cells(cell,1)-faces(face,1))^2*(cells(cell,2)-faces(face,2))^3;
        D{face}(j,21)=(cells(cell,1)-faces(face,1))^1*(cells(cell,2)-faces(face,2))^4;
        
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
            facestencil=stencil_faces(face,j);           
            if ponderado==true
                if facestencil~=face
                    dx=(faces(facestencil,1)-faces(face,1));
                    dy=(faces(facestencil,2)-faces(face,2));
                else
                    guarda(face)=n_cell+j;
                    aux=face_cells(face,1);
                    dx=(cells(aux,1)-faces(face,1));
                    dy=(cells(aux,2)-faces(face,2));
%                      dx=dx/10;
%                      dy=dy/10;
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
            
      
      
  if  face_bound(face,4)==1        
           for q=1:3

            Q_border{face}(q,1)=0;
            %
           Q_border{face}(q,2)=normal_x;
            Q_border{face}(q,3)=normal_y;
            %
           Q_border{face}(q,4)=2*(x(q)-faces(face,1))*normal_x;
            Q_border{face}(q,5)=2*(y(q)-faces(face,2))*normal_y;
            Q_border{face}(q,6)=(x(q)-faces(face,1))*normal_y+(y(q)-faces(face,2))*normal_x;
            %
            Q_border{face}(q,7)=(3*(x(q)-faces(face,1))^2)*normal_x;
            Q_border{face}(q,8)=(3*(y(q)-faces(face,2))^2)*normal_y;
            Q_border{face}(q,9)=2*(x(q)-faces(face,1))*(y(q)-faces(face,2))*normal_x+((x(q)-faces(face,1))^2)*normal_y;
            Q_border{face}(q,10)=((y(q)-faces(face,2))^2)*normal_x+2*(x(q)-faces(face,1))*(y(q)-faces(face,2))*normal_y;
            %
            Q_border{face}(q,11)=4*(x(q)-faces(face,1))^3*normal_x;
            Q_border{face}(q,12)=4*(y(q)-faces(face,2))^3*normal_y;           
            Q_border{face}(q,13)=3*(x(q)-faces(face,1))^2*(y(q)-faces(face,2))*normal_x+(x(q)-faces(face,1))^3*normal_y;
            Q_border{face}(q,14)=2*(x(q)-faces(face,1))*(y(q)-faces(face,2))^2*normal_x+2*(x(q)-faces(face,1))^2*(y(q)-faces(face,2))*normal_y;
            Q_border{face}(q,15)=(y(q)-faces(face,2))^3*normal_x+3*(x(q)-faces(face,1))^1*(y(q)-faces(face,2))^2*normal_y;

            %
            Q_border{face}(q,16)=5*(x(q)-faces(face,1))^4*normal_x;
            Q_border{face}(q,17)=5*(y(q)-faces(face,2))^4*normal_y;
            Q_border{face}(q,18)=4*(x(q)-faces(face,1))^3*(y(q)-faces(face,2))^1*normal_x+(x(q)-faces(face,1))^4*normal_y;
            Q_border{face}(q,19)=3*(x(q)-faces(face,1))^2*(y(q)-faces(face,2))^2*normal_x+2*(x(q)-faces(face,1))^3*(y(q)-faces(face,2))^1*normal_y;
            Q_border{face}(q,20)=2*(x(q)-faces(face,1))^1*(y(q)-faces(face,2))^3*normal_x+3*(x(q)-faces(face,1))^2*(y(q)-faces(face,2))^2*normal_y;
            Q_border{face}(q,21)=(y(q)-faces(face,2))^4*normal_x+4*(x(q)-faces(face,1))^1*(y(q)-faces(face,2))^3*normal_y;
            

            %

    
            c{face}(q,1)=SolutionDiffusion(solution,x(q),y(q),'xflux')*normal_x;
    end
  end
  
  
  %% robin
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
            
            
         
%%%%%%%%%%
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
            
            
         if  face_bound(face,4)==2        
           for q=1:3

             Q_border{face}(q,1)=0;
            %
             Q_border{face}(q,2)=1*normal_x;
             Q_border{face}(q,3)=0;
            %
             Q_border{face}(q,4)=2*(x(q)-faces(face,1))*normal_x;
             Q_border{face}(q,5)=0;
             Q_border{face}(q,6)=(y(q)-faces(face,2))*normal_x;
            %
            Q_border{face}(q,7)=3*(x(q)-faces(face,1))^2*normal_x;
             Q_border{face}(q,8)=0;
            Q_border{face}(q,9)=2*(x(q)-faces(face,1))^1*(y(q)-faces(face,2))*normal_x;
             Q_border{face}(q,10)=(y(q)-faces(face,2))^2*normal_x;
            %
             Q_border{face}(q,11)=4*(x(q)-faces(face,1))^3*normal_x;
            Q_border{face}(q,12)=4*(y(q)-faces(face,2))^3*normal_y;           
            Q_border{face}(q,13)=3*(x(q)-faces(face,1))^2*(y(q)-faces(face,2))*normal_x+(x(q)-faces(face,1))^3*normal_y;
            Q_border{face}(q,14)=2*(x(q)-faces(face,1))*(y(q)-faces(face,2))^2*normal_x+2*(x(q)-faces(face,1))^2*(y(q)-faces(face,2))*normal_y;
            Q_border{face}(q,15)=(y(q)-faces(face,2))^3*normal_x+3*(x(q)-faces(face,1))^1*(y(q)-faces(face,2))^2*normal_y;

            %
            Q_border{face}(q,16)=5*(x(q)-faces(face,1))^4*normal_x;
            Q_border{face}(q,17)=5*(y(q)-faces(face,2))^4*normal_y;
            Q_border{face}(q,18)=4*(x(q)-faces(face,1))^3*(y(q)-faces(face,2))^1*normal_x+(x(q)-faces(face,1))^4*normal_y;
            Q_border{face}(q,19)=3*(x(q)-faces(face,1))^2*(y(q)-faces(face,2))^2*normal_x+2*(x(q)-faces(face,1))^3*(y(q)-faces(face,2))^1*normal_y;
            Q_border{face}(q,20)=2*(x(q)-faces(face,1))^1*(y(q)-faces(face,2))^3*normal_x+3*(x(q)-faces(face,1))^2*(y(q)-faces(face,2))^2*normal_y;
            Q_border{face}(q,21)=(y(q)-faces(face,2))^4*normal_x+4*(x(q)-faces(face,1))^1*(y(q)-faces(face,2))^3*normal_y;
            

            
       
            
            Q_border{face}(q,1)=(gamma_diff/u_convec_x)*Q_border{face}(q,1)+area_robin*1;
            %
             Q_border{face}(q,2)=(gamma_diff/u_convec_x)*Q_border{face}(q,2)+area_robin*(x(q)-faces(face,1));
             Q_border{face}(q,3)=(gamma_diff/u_convec_x)*Q_border{face}(q,3)+area_robin*(y(q)-faces(face,2));
            %
             Q_border{face}(q,4)=(gamma_diff/u_convec_x)*Q_border{face}(q,4)+area_robin*(x(q)-faces(face,1))^2;
             Q_border{face}(q,5)=(gamma_diff/u_convec_x)*Q_border{face}(q,5)+area_robin*(y(q)-faces(face,2))^2;
             Q_border{face}(q,6)=(gamma_diff/u_convec_x)*Q_border{face}(q,6)+area_robin*(x(q)-faces(face,1))*(y(q)-faces(face,2));
            %
            Q_border{face}(q,7)=(gamma_diff/u_convec_x)*Q_border{face}(q,7)+area_robin*(x(q)-faces(face,1))^3;
             Q_border{face}(q,8)=(gamma_diff/u_convec_x)*Q_border{face}(q,8)+area_robin*(y(q)-faces(face,2))^3;
            Q_border{face}(q,9)=(gamma_diff/u_convec_x)*Q_border{face}(q,9)+area_robin*(x(q)-faces(face,1))^2*(y(q)-faces(face,2));
             Q_border{face}(q,10)=(gamma_diff/u_convec_x)*Q_border{face}(q,10)+area_robin*(x(q)-faces(face,1))*(y(q)-faces(face,2))^2;
    
              Q_border{face}(q,11)=(gamma_diff/u_convec_x)*Q_border{face}(q,11)+area_robin*(x(q)-faces(face,1))^4;
            Q_border{face}(q,12)=(gamma_diff/u_convec_x)*Q_border{face}(q,12)+area_robin*(y(q)-faces(face,2))^4;
             Q_border{face}(q,13)=(gamma_diff/u_convec_x)*Q_border{face}(q,13)+area_robin*(x(q)-faces(face,1))^3*(y(q)-faces(face,2))^1;
             Q_border{face}(q,14)=(gamma_diff/u_convec_x)*Q_border{face}(q,14)+area_robin*(x(q)-faces(face,1))^2*(y(q)-faces(face,2))^2;
             Q_border{face}(q,15)=(gamma_diff/u_convec_x)*Q_border{face}(q,15)+area_robin*(x(q)-faces(face,1))^1*(y(q)-faces(face,2))^3;
            %
            Q_border{face}(q,16)=(gamma_diff/u_convec_x)*Q_border{face}(q,16)+area_robin*(x(q)-faces(face,1))^5;
           Q_border{face}(q,17)=(gamma_diff/u_convec_x)*Q_border{face}(q,17)+area_robin*(y(q)-faces(face,2))^5;
            Q_border{face}(q,18)=(gamma_diff/u_convec_x)*Q_border{face}(q,18)+area_robin*(x(q)-faces(face,1))^4*(y(q)-faces(face,2))^1;
             Q_border{face}(q,19)=(gamma_diff/u_convec_x)*Q_border{face}(q,19)+area_robin*(x(q)-faces(face,1))^3*(y(q)-faces(face,2))^2;
             Q_border{face}(q,20)=(gamma_diff/u_convec_x)*Q_border{face}(q,20)+area_robin*(x(q)-faces(face,1))^2*(y(q)-faces(face,2))^3;
             Q_border{face}(q,21)=(gamma_diff/u_convec_x)*Q_border{face}(q,21)+area_robin*(x(q)-faces(face,1))^1*(y(q)-faces(face,2))^4;
             
             
             
            c{face}(q,1)=(gamma_diff/u_convec_x)*SolutionDiffusion(solution,x(q),y(q),'xflux')*normal_x+area_robin*SolutionDiffusion(solution,x(q),y(q),'anal');
    end
   end 
  
  
  
  
  
  
  
   
  %% dirichlet
  
  else
            soluti{face}(n_cell+j) = SolutionDiffusion(solution,faces(facestencil,1),faces(facestencil,2),'anal'); 
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
            
    if  face_bound(face,4)==0            
    for q=1:3

             Q_border{face}(q,1)=1;
            %
             Q_border{face}(q,2)=(x(q)-faces(face,1));
             Q_border{face}(q,3)=(y(q)-faces(face,2));
            %
             Q_border{face}(q,4)=(x(q)-faces(face,1))^2;
             Q_border{face}(q,5)=(y(q)-faces(face,2))^2;
             Q_border{face}(q,6)=(x(q)-faces(face,1))*(y(q)-faces(face,2));
            %
            Q_border{face}(q,7)=(x(q)-faces(face,1))^3;
             Q_border{face}(q,8)=(y(q)-faces(face,2))^3;
            Q_border{face}(q,9)=(x(q)-faces(face,1))^2*(y(q)-faces(face,2));
             Q_border{face}(q,10)=(x(q)-faces(face,1))*(y(q)-faces(face,2))^2;
            %
             Q_border{face}(q,11)=(x(q)-faces(face,1))^4;
            Q_border{face}(q,12)=(y(q)-faces(face,2))^4;
             Q_border{face}(q,13)=(x(q)-faces(face,1))^3*(y(q)-faces(face,2))^1;
             Q_border{face}(q,14)=(x(q)-faces(face,1))^2*(y(q)-faces(face,2))^2;
             Q_border{face}(q,15)=(x(q)-faces(face,1))^1*(y(q)-faces(face,2))^3;
            %
            Q_border{face}(q,16)=(x(q)-faces(face,1))^5;
           Q_border{face}(q,17)=(y(q)-faces(face,2))^5;
            Q_border{face}(q,18)=(x(q)-faces(face,1))^4*(y(q)-faces(face,2))^1;
             Q_border{face}(q,19)=(x(q)-faces(face,1))^3*(y(q)-faces(face,2))^2;
             Q_border{face}(q,20)=(x(q)-faces(face,1))^2*(y(q)-faces(face,2))^3;
             Q_border{face}(q,21)=(x(q)-faces(face,1))^1*(y(q)-faces(face,2))^4;
             
             
              if extended_stencil && nx==0       
        Q_border{face}(q,22)=(y(q)-faces(face,2))^6;
        Q_border{face}(q,23)=(x(q)-faces(face,1))^5*(y(q)-faces(face,2))^1;
        Q_border{face}(q,24)=(x(q)-faces(face,1))^4*(y(q)-faces(face,2))^2;
       Q_border{face}(q,25)=(x(q)-faces(face,1))^3*(y(q)-faces(face,2))^3;
        Q_border{face}(q,26)=(x(q)-faces(face,1))^2*(y(q)-faces(face,2))^4;
        Q_border{face}(q,27)=(x(q)-faces(face,1))^1*(y(q)-faces(face,2))^5;
        elseif extended_stencil && ny==0
        Q_border{face}(q,22)=(x(q)-faces(face,1))^6;
       Q_border{face}(q,23)=(x(q)-faces(face,1))^5*(y(q)-faces(face,2))^1;
        Q_border{face}(q,24)=(x(q)-faces(face,1))^4*(y(q)-faces(face,2))^2;
        Q_border{face}(q,25)=(x(q)-faces(face,1))^3*(y(q)-faces(face,2))^3;
        Q_border{face}(q,26)=(x(q)-faces(face,1))^2*(y(q)-faces(face,2))^4;
        Q_border{face}(q,27)=(x(q)-faces(face,1))^1*(y(q)-faces(face,2))^5;
        end
    
            c{face}(q,1)=SolutionDiffusion(solution,x(q),y(q),'anal');
    end
    end
    
    %% Determinação da Matriz T %%
    % D*c=phi -> c=inv(D'D)*D'*phi -> T=inv(D'D)*D' %

  if  face_bound(face,4)==0 ||   face_bound(face,4)==1 ||   face_bound(face,4)==2          
    D{face}(guarda(face),:)=[];
    D1{face}(guarda(face),:)=[]; 
    
    Q_border{face}=Q_border{face}';
   
    
    T1{face} = inv(D{face}'*D1{face});
    
    M{face} = inv(D{face}'*D1{face})*D1{face}';
   
    K{face} = T1{face}*Q_border{face}*inv(Q_border{face}'*T1{face}*Q_border{face});
    
    
    T_border{face}= M{face}- K{face}*Q_border{face}'*M{face};
    
    constrained_source{face}= K{face}*c{face};
    
    
%     
%     restos{face}{q}= sum(sum(inv(D_border{face}{q}'*D1_border{face}{q})*(D_border{face}{q}'*D1_border{face}{q})-eye(21)));
%     
  end
    end
    end
end
end

% testar_constraints(M,T_border,constrained_source,c,Q_border,soluti,6);

 pause(0.1)
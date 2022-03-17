function [T,D]=Reconstruction2ndOrder(stencil_cells,stencil_faces,stencil_size,ponderado)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   26 de Setembro de 2016                                          %
%                                                                                                   %
% Função que implementa a Solução Numérica a partir do Minimos Quadrados de 2ª Ordem                %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global L Lref cell_side vert_side cell_num face_num vert_num w;
global verts cells cell_verts cell_faces cell_vol faces face_area face_vert cell_norm;
global face_bound cell_bound face_cells vert_cells u_convec_x u_convec_y gamma_diff;
global phi lap_phi  phi_faces flux_phi_faces extended_stencil neuman dimensional_correction robin;
%
for i=1:face_num
    %%% _> fft
    
    face=i;
    
%     % > debug.
%     hold on;
%     plot(verts(:,1),verts(:,2),'-ok');
    
    % Normal Exterior da Face %
    vert=face_vert(face,:);
    nx=verts(vert(1),1)-verts(vert(2),1);
    ny=verts(vert(1),2)-verts(vert(2),2);
    
%     % > debug.
%     hold on;
%     plot(verts(vert(1),1),verts(vert(1),2),'-ob');
%     plot(verts(vert(2),1),verts(vert(2),2),'-ob');
    
    if nx~=0
        nx=1;
    elseif ny~=0
        ny=1;
    end
    %%% _> fft
    
    
    
    
    
    %% Construção da Matriz D para cada Face, Celulas do Stencil %%
    n_cell=stencil_size(face,2);
    
    for j=1:n_cell
        celula=stencil_cells(face,j);
        
        if ponderado
            dx=(cells(celula,1)-faces(face,1));
            dy=(cells(celula,2)-faces(face,2));
            w=1/(sqrt(dx^2+dy^2))^2;
        else
            w=1;
        end
 
        % Matriz D tendo em conta a ponderação %
        D1{face}(j,1)=w*1;
        D1{face}(j,2)=w*(cells(celula,1)-faces(face,1));
        D1{face}(j,3)=w*(cells(celula,2)-faces(face,2));

        if extended_stencil && nx==0
            D1{face}(j,4)=w*(cells(celula,2)-faces(face,2))^2;
            D1{face}(j,5)=w*(cells(celula,1)-faces(face,1))*(cells(celula,2)-faces(face,2));
        elseif extended_stencil && ny==0
            D1{face}(j,4)=w*(cells(celula,1)-faces(face,1))^2;
            D1{face}(j,5)=w*(cells(celula,1)-faces(face,1))*(cells(celula,2)-faces(face,2));
        end
                
        % Matriz D %
        D{face}(j,1)=1;
        D{face}(j,2)=(cells(celula,1)-faces(face,1));
        D{face}(j,3)=(cells(celula,2)-faces(face,2));
        
        if extended_stencil && nx==0
            D{face}(j,4)=(cells(celula,2)-faces(face,2))^2;
            D{face}(j,5)=(cells(celula,1)-faces(face,1))*(cells(celula,2)-faces(face,2));
        elseif extended_stencil && ny==0
            D{face}(j,4)=(cells(celula,1)-faces(face,1))^2;
            D{face}(j,5)=(cells(celula,1)-faces(face,1))*(cells(celula,2)-faces(face,2));
        end
    end
    %
    %% Construção da Matriz D para cada Face, Faces do Stencil %%
    n_face=stencil_size(face,3);
    if n_face~=0
        for j=1:n_face
            facestencil=stencil_faces(face,j);
            %
            if ponderado
                if facestencil~=face
                    dx=(faces(facestencil,1)-faces(face,1));
                    dy=(faces(facestencil,2)-faces(face,2));
                else
                    aux=face_cells(face,1);
                    dx=(cells(aux,1)-faces(face,1));
                    dy=(cells(aux,2)-faces(face,2));
                    %                     dx=0.000000001;
                    %                     dy=0.000000001;
                end
                w=1/(sqrt(dx^2+dy^2))^2;
            else
                w=1;
            end

            %%  NEUMAN
            
            normal_x=face_bound(facestencil,2)*face_area(facestencil,1);
            normal_y=face_bound(facestencil,3)*face_area(facestencil,1);
            area_robin=1*face_area(facestencil,1);
            %
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
                D1{face}(n_cell+j,2)=w*normal_x;
                D1{face}(n_cell+j,3)=0;

                if extended_stencil && nx==0
                    D1{face}(n_cell+j,4)=2*w*(faces(facestencil,2)-faces(face,2))*normal_y;
                    D1{face}(n_cell+j,5)=w*(faces(facestencil,2)-faces(face,2))*normal_x+w*(faces(facestencil,1)-faces(face,1))*normal_y;
                elseif extended_stencil && ny==0
                    D1{face}(n_cell+j,4)=2*w*(faces(facestencil,1)-faces(face,1))*normal_x;
                    D1{face}(n_cell+j,5)=w*(faces(facestencil,2)-faces(face,2))*normal_x+w*(faces(facestencil,1)-faces(face,1))*normal_y;
                end
                
                % Matriz D %
                D{face}(n_cell+j,1)=0;
                D{face}(n_cell+j,2)=normal_x;
                D{face}(n_cell+j,3)=0;                
                
                if extended_stencil && nx==0
                    D{face}(n_cell+j,4)=2*(faces(facestencil,2)-faces(face,2))*normal_y;
                    D{face}(n_cell+j,5)=(faces(facestencil,2)-faces(face,2))*normal_x+(faces(facestencil,1)-faces(face,1))*normal_y;
                elseif extended_stencil && ny==0
                    D{face}(n_cell+j,4)=2*(faces(facestencil,1)-faces(face,1))*normal_x;
                    D{face}(n_cell+j,5)=(faces(facestencil,2)-faces(face,2))*normal_x+(faces(facestencil,1)-faces(face,1))*normal_y;
                end
                
                %% ROBIN                
            elseif robin && face_bound(facestencil,4)==2
                
                D1{face}(n_cell+j,1)=0;
                
                D1{face}(n_cell+j,2)=w*normal_x;
                D1{face}(n_cell+j,3)=0;

                % Matriz D %
                D{face}(n_cell+j,1)=0;
                D{face}(n_cell+j,2)=1*normal_x;
                D{face}(n_cell+j,3)=0;

                % Matriz D tendo em conta a ponderação %
                D1{face}(n_cell+j,1)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,1)+area_robin*w*1;
                %
                D1{face}(n_cell+j,2)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,2)+area_robin*w*(faces(facestencil,1)-faces(face,1));
                D1{face}(n_cell+j,3)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,3)+area_robin*w*(faces(facestencil,2)-faces(face,2));

                % Matriz D %
                D{face}(n_cell+j,1)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,1)+area_robin*1;
                %
                D{face}(n_cell+j,2)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,2)+area_robin*(faces(facestencil,1)-faces(face,1));
                D{face}(n_cell+j,3)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,3)+area_robin*(faces(facestencil,2)-faces(face,2));
                
 
                %% Drichlet
            else
                % Matriz D tendo em conta a ponderação %
                D1{face}(n_cell+j,1)=w*1;
                %
                D1{face}(n_cell+j,2)=w*(faces(facestencil,1)-faces(face,1));
                D1{face}(n_cell+j,3)=w*(faces(facestencil,2)-faces(face,2));
                %
                
                if extended_stencil && nx==0
                    D1{face}(n_cell+j,4)=w*(faces(facestencil,2)-faces(face,2))^2;
                    D1{face}(n_cell+j,5)=w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2));
                elseif extended_stencil && ny==0
                    D1{face}(n_cell+j,4)=w*(faces(facestencil,1)-faces(face,1))^2;
                    D1{face}(n_cell+j,5)=w*(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2));
                end
                
                % Matriz D %
                D{face}(n_cell+j,1)=1;
                %
                D{face}(n_cell+j,2)=(faces(facestencil,1)-faces(face,1));
                D{face}(n_cell+j,3)=(faces(facestencil,2)-faces(face,2));
                                
                if extended_stencil && nx==0
                    D{face}(n_cell+j,4)=(faces(facestencil,2)-faces(face,2))^2;
                    D{face}(n_cell+j,5)=(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2));
                elseif extended_stencil && ny==0
                    D{face}(n_cell+j,4)=(faces(facestencil,1)-faces(face,1))^2;
                    D{face}(n_cell+j,5)=(faces(facestencil,1)-faces(face,1))*(faces(facestencil,2)-faces(face,2));
                end
            end
        end
    end
    %
    %% Determinação da Matriz T %%
    %
    % D*c=phi -> c=inv(D'D)*D'*phi -> T=inv(D'D)*D' %
    T{face}=inv(D{face}'*D1{face})*D1{face}';
end
pause(0.0001)
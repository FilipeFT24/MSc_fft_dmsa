function [T_border,D, constrained_source]=Reconstruction2ndOrder_borders_v2(stencil_cells,stencil_faces,stencil_size,ponderado)
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
global face_bound cell_bound face_cells vert_cells robin u_convec_x u_convec_y gamma_diff;
global phi lap_phi  phi_faces flux_phi_faces solution G extended_stencil face_w_gauss dimensional_correction neuman;
%
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
    %% Construção da Matriz D para cada Face, Celulas do Stencil %%
    n_cell=stencil_size(face,2);
    %
    for j=1:n_cell
        celula=stencil_cells(face,j);
        %
        if ponderado
            dx=(cells(celula,1)-faces(face,1));
            dy=(cells(celula,2)-faces(face,2));
            w=1/(sqrt(dx^2+dy^2))^2;
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
          if extended_stencil && nx==0       
        D1{face}(j,4)=w*(cells(celula,2)-faces(face,2))^2;
        D1{face}(j,5)=w*(cells(celula,1)-faces(face,1))*(cells(celula,2)-faces(face,2));
        elseif extended_stencil && ny==0
        D1{face}(j,4)=w*(cells(celula,1)-faces(face,1))^2;
        D1{face}(j,5)=w*(cells(celula,1)-faces(face,1))*(cells(celula,2)-faces(face,2));
        end
        
        % Matriz D %
        D{face}(j,1)=1;
        %
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
                    guarda(face)=n_cell+j;
                    aux=face_cells(face,1);
                    dx=(cells(aux,1)-faces(face,1));
                    dy=(cells(aux,2)-faces(face,2));
%                     dx=1000000;
%                     dy=1000000;
                end
                w=1/(sqrt(dx^2+dy^2))^2;
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
        D1{face}(n_cell+j,3)=0;
        %
  
             % Matriz D %
        D{face}(n_cell+j,1)=0;
        %
        D{face}(n_cell+j,2)=1*normal_x;
        D{face}(n_cell+j,3)=0;


   if  face_bound(face,4)==1        
           for q=1:1

             Q_border{face}(q,1)=0;
            %
             Q_border{face}(q,2)=1*normal_x;
             Q_border{face}(q,3)=0;
            %

            %

    
            c{face}(q,1)=SolutionDiffusion(solution,x(q),y(q),'xflux')*normal_x;
    end
   end       
  %% robin              
 elseif robin && face_bound(facestencil,4)==2
                
%             if ~dimensional_correction
%               normal_x=normal_x/face_area(facestencil,1);
%               normal_y=normal_y/face_area(facestencil,1);
%              end    
                
                
                
            D1{face}(n_cell+j,1)=0;
            %
            D1{face}(n_cell+j,2)=w*normal_x;
            D1{face}(n_cell+j,3)=0;
            %
            
      
            
            % Matriz D %
            D{face}(n_cell+j,1)=0;
            %
            D{face}(n_cell+j,2)=normal_x;
            D{face}(n_cell+j,3)=0;
            
            
      
        
        
        
        % Matriz D tendo em conta a ponderação %
            D1{face}(n_cell+j,1)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,1)+area_robin*w*1;
            %
            D1{face}(n_cell+j,2)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,2)+area_robin*w*(faces(facestencil,1)-faces(face,1));
            D1{face}(n_cell+j,3)=(gamma_diff/u_convec_x)*D1{face}(n_cell+j,3)+area_robin*w*(faces(facestencil,2)-faces(face,2));
            %
            
        
            
            % Matriz D %
            D{face}(n_cell+j,1)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,1)+area_robin*1;
            %
            D{face}(n_cell+j,2)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,2)+area_robin*(faces(facestencil,1)-faces(face,1));
            D{face}(n_cell+j,3)=(gamma_diff/u_convec_x)*D{face}(n_cell+j,3)+area_robin*(faces(facestencil,2)-faces(face,2));
            
            
      if  face_bound(face,4)==2        
           for q=1:1

             Q_border{face}(q,1)=0;
            %
             Q_border{face}(q,2)=1*normal_x;
             Q_border{face}(q,3)=0;
            %
             
            %
       
            
            Q_border{face}(q,1)=(gamma_diff/u_convec_x)*Q_border{face}(q,1)+area_robin*1;
            %
             Q_border{face}(q,2)=(gamma_diff/u_convec_x)*Q_border{face}(q,2)+area_robin*(x(q)-faces(face,1));
             Q_border{face}(q,3)=(gamma_diff/u_convec_x)*Q_border{face}(q,3)+area_robin*(y(q)-faces(face,2));
            %
            
            c{face}(q,1)=(gamma_diff/u_convec_x)*SolutionDiffusion(solution,x(q),y(q),'xflux')*normal_x+area_robin*SolutionDiffusion(solution,x(q),y(q),'anal');
    end
   end  
        
           
   
  
      
           
 %% dirichlet          
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

        
  if  face_bound(face,4)==0       
     for q=1:1

             Q_border{face}(q,1)=1;
            %
             Q_border{face}(q,2)=(x(q)-faces(face,1));
             Q_border{face}(q,3)=(y(q)-faces(face,2));
            %    
            
        if extended_stencil && nx==0       
        Q_border{face}(q,4)=(y(q)-faces(face,2))^2;
       Q_border{face}(q,5)=(x(q)-faces(face,1))*(y(q)-faces(face,2));
        elseif extended_stencil && ny==0
        Q_border{face}(q,4)=(x(q)-faces(face,1))^2;
        Q_border{face}(q,5)=(x(q)-faces(face,1))*(y(q)-faces(face,2));
        end
            
            c{face}(q,1)=SolutionDiffusion(solution,x(q),y(q),'anal');
     end
  end       
        end
        

    %
    %% Determinação da Matriz T %%
    %
    % D*c=phi -> c=inv(D'D)*D'*phi -> T=inv(D'D)*D' %
    if  face_bound(face,4)==0 ||   face_bound(face,4)==1  ||   face_bound(face,4)==2         
    D{face}(guarda(face),:)=[];
    D1{face}(guarda(face),:)=[]; 
    Q_border{face}=Q_border{face}';
    
    T1{face} = inv(D{face}'*D1{face});
    
    M{face} = inv(D{face}'*D1{face})*D1{face}';
   
    K{face} = T1{face}*Q_border{face}*inv(Q_border{face}'*T1{face}*Q_border{face});    
    
    T_border{face}= M{face}- K{face}*Q_border{face}'*M{face};
    
    constrained_source{face}= K{face}*c{face};
    end
    end
      end
      end
end

pause(0.1)
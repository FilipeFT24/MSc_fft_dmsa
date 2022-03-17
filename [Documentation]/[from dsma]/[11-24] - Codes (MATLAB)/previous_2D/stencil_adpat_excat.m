function []=stencil_adpat_excat()

global TempoGlobal tempo_malha tempo_anal tempo_num tempo_A tempo_stencil tempo_rec tempo_gmres tempo_erro tempo_plots;
global fid;
global L Lref cell_side vert_side cell_num face_num vert_num;
global verts cells faces;
global cell_verts cell_faces cell_vol cell_norm cell_bound cell_vert_num cell_face_num;
global face_vert face_cells face_area face_bound;
global vert_cells vert_cell_num vert_face_num;
global phi lap_phi  phi_faces flux_phi_faces;
global phi_num lap_phi_num A source source_faces source_cells stencil_cells stencil_faces stencil_size T D;
global norma1_phi erro_phi_max erro_phi norma1_lap erro_lap_max erro_lap X;
global u_convec_x u_convec_y gamma_diff order fig face_w_gauss mix_method G solution T_border restos constrained_source extended_stencil increase_gauss_points;



stencil_cells_adapt=stencil_cells;
stencil_faces_adapt=stencil_faces;
stencil_size_adapt= stencil_size;


for face_v2 = 1:face_num
    
  %% caso dos cantos  
%    if stencil_size(face_v2,1)< order*(order+1)
       if ismember(1,stencil_cells(face_v2,:))
           stencil_faces_adapt(face_v2,stencil_size_adapt(face_v2,3)+1)= face_num+1; %canto inf. esq.
           stencil_size_adapt(face_v2,1)=stencil_size_adapt(face_v2,1)+1;
           stencil_size_adapt(face_v2,3)=stencil_size_adapt(face_v2,3)+1;
       end
       
       if ismember(cell_side,stencil_cells(face_v2,:))
           stencil_faces_adapt(face_v2,stencil_size_adapt(face_v2,3)+1)= face_num+2; %canto sup. esq.
           stencil_size_adapt(face_v2,1)=stencil_size_adapt(face_v2,1)+1;
           stencil_size_adapt(face_v2,3)=stencil_size_adapt(face_v2,3)+1;
       end
       
        if ismember(cell_num+1-cell_side,stencil_cells(face_v2,:))
            stencil_faces_adapt(face_v2,stencil_size_adapt(face_v2,3)+1)= face_num+3; %canto inf. dto.
           stencil_size_adapt(face_v2,1)=stencil_size_adapt(face_v2,1)+1;
           stencil_size_adapt(face_v2,3)=stencil_size_adapt(face_v2,3)+1;
        end
       
         if ismember(cell_num,stencil_cells(face_v2,:))
           stencil_faces_adapt(face_v2,stencil_size_adapt(face_v2,3)+1)= face_num+4;%canto sup. dto.
           stencil_size_adapt(face_v2,1)=stencil_size_adapt(face_v2,1)+1;
           stencil_size_adapt(face_v2,3)=stencil_size_adapt(face_v2,3)+1;
       end
       
%    end
   
%% numero de pontos em x vs y
if stencil_size_adapt(face_v2,1)>order*(order+1)
 aux_y=find((stencil_faces_adapt(face_v2,:)<= cell_side*(cell_side+1) | stencil_faces_adapt(face_v2,:)>face_num) & stencil_faces_adapt(face_v2,:)~=0);
 aux_x=find((stencil_faces_adapt(face_v2,:)> cell_side*(cell_side+1) | stencil_faces_adapt(face_v2,:)>face_num) & stencil_faces_adapt(face_v2,:)~=0);

 tamanho_x =size(aux_x,2);
 tamanho_y =size(aux_y,2);

  vert=face_vert(face_v2,:);
  nx=verts(vert(1),1)-verts(vert(2),1);
  ny=verts(vert(1),2)-verts(vert(2),2); 
  
  
  alert_x=0;
  alert_y=0;
  
  if (nx==0 && tamanho_x ~=order) || (nx~=0 && tamanho_x~=order+1)
    alert_x=1;  
  end
  if (ny==0 && tamanho_y ~=order) || (ny~=0 && tamanho_y~=order+1)
     alert_y=1; 
  end

   %% caso de excesso de faces fronteira  
    
   
   
   if alert_x==1 
      aux1= find((stencil_faces_adapt(face_v2,:)<= cell_side*(cell_side+1) | stencil_faces_adapt(face_v2,:)>face_num) & stencil_faces_adapt(face_v2,:)~=0);
      stencil_faces_adapt_copy= stencil_faces_adapt(face_v2,:);
      stencil_faces_adapt_copy(aux1)=[];
      stencil_faces_adapt_copy(numel(stencil_faces_adapt(face_v2,:))) = 0;
      stencil_faces_adapt(face_v2,:)=stencil_faces_adapt_copy;
      stencil_size_adapt(face_v2,1)=stencil_size_adapt(face_v2,1)-size(aux1,2);
      stencil_size_adapt(face_v2,3)=stencil_size_adapt(face_v2,3)-size(aux1,2);     
   end
   
   
   
   if alert_y==1
%        if face_v2==245
%     pause(0.1)
% end
       aux1= find((stencil_faces_adapt(face_v2,:)> cell_side*(cell_side+1) | stencil_faces_adapt(face_v2,:)>face_num) & stencil_faces_adapt(face_v2,:)~=0);
      stencil_faces_adapt_copy= stencil_faces_adapt(face_v2,:);
      stencil_faces_adapt_copy(aux1)=[];
      stencil_faces_adapt_copy(numel(stencil_faces_adapt(face_v2,:))) = 0;
      stencil_faces_adapt(face_v2,:)=stencil_faces_adapt_copy;
      stencil_size_adapt(face_v2,1)=stencil_size_adapt(face_v2,1)-size(aux1,2);
      stencil_size_adapt(face_v2,3)=stencil_size_adapt(face_v2,3)-size(aux1,2);
   end
      

% 
% %% casos de excesso de faces fronteira 2
% 
% if stencil_size_adapt(face_v2,1)> order*(order+1)
%       aux1= find((stencil_faces_adapt(face_v2,:)~=0));
%       stencil_faces_adapt_copy= stencil_faces_adapt(face_v2,:);
%       stencil_faces_adapt_copy(aux1)=0;
%       stencil_faces_adapt(face_v2,:)=stencil_faces_adapt_copy;
%       stencil_size_adapt(face_v2,1)=stencil_size_adapt(face_v2,1)-stencil_size_adapt(face_v2,3);
%       stencil_size_adapt(face_v2,3)=stencil_size_adapt(face_v2,3)-stencil_size_adapt(face_v2,3);  
% end
end
end
%% Adicionar pontos do cantos ao vetor faces

faces(face_num+1,1)=0;
faces(face_num+1,2)=0;
phi_faces(face_num+1)= SolutionDiffusion(solution,0,0,'anal');

faces(face_num+2,1)=0;
faces(face_num+2,2)=L;
phi_faces(face_num+2)= SolutionDiffusion(solution,0,L,'anal');

faces(face_num+3,1)=L;
faces(face_num+3,2)=0;
phi_faces(face_num+3)= SolutionDiffusion(solution,L,0,'anal');

faces(face_num+4,1)=L;
faces(face_num+4,2)=L;
phi_faces(face_num+4)= SolutionDiffusion(solution,L,L,'anal');


%% Alterar o vetor stencil

stencil_cells=stencil_cells_adapt;
stencil_faces=stencil_faces_adapt;
 stencil_size=stencil_size_adapt;

pause(0.1)
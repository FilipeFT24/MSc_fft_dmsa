function []=stencil_visualization()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       Filipe Diogo                                                %
%                                         07/03/2018                                               %
%                                                                                                   %
%                                                                                                   %
                                                                                                      %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



global L Lref cell_side vert_side cell_num face_num vert_num w;
global verts cells cell_verts cell_faces cell_vol faces face_area face_vert cell_norm;
global face_bound cell_bound face_cells vert_cells stencil_size stencil_cells stencil_faces T;
global phi lap_phi  phi_faces flux_phi_faces phi_num order fig face_w_gauss neuman T_border G constrained_source extended_stencil;


%  close(fig)
figure()

for face_v2 = 1:face_num
    
clear cell_aux_v2 solucao_analitica solucao_numerica localizacao_celulas_stencil C_v2 polinomio phi_estimado peso volume_celula volume_total
    
xc=faces(face_v2,1);
yc=faces(face_v2,2);

guarda_face=0;
localizacao_celulas_stencil_neuman=[];

% aux_yc(face_v2)= yc;
% aux_xc(face_v2)=xc;

 for k=1:stencil_size(face_v2,2)
                %
                % Celulas do Stencil %
                 cell_aux_v2(k)=stencil_cells(face_v2,k);  
%                 solucao_analitica(k) = phi(1,cell_aux_v2(k));
%                 solucao_numerica(k) = phi_num(cell_aux_v2(k),1);
                localizacao_celulas_stencil(k,:) = [cells(cell_aux_v2(k),1) cells(cell_aux_v2(k),2)];  
 end    
 
  for cont=1:stencil_size(face_v2,3)
                %
                % faces do Stencil %
%                 if stencil_faces(face_v2,cont) ~= face_v2
                 cell_aux_v2(k+cont)=stencil_faces(face_v2,cont);
%                 solucao_analitica(k+cont) = phi_faces(cell_aux_v2(k+cont));
%                 solucao_numerica(k+cont) = phi_faces(cell_aux_v2(k+cont));
                localizacao_celulas_stencil(k+cont,:) = [faces(cell_aux_v2(k+cont),1) faces(cell_aux_v2(k+cont),2)];
if stencil_faces(face_v2,cont) == face_v2
 guarda_face = k+cont;
%   cell_aux_v2(k+cont)=face_cells(face_v2,1);
end
%                 else
%                 guarda_face = k+cont;
%                 cell_aux_v2(k+cont)=face_cells(face_v2,1); 
%                 solucao_analitica(k+cont) = phi(1,cell_aux_v2(k+cont));
%                 solucao_numerica(k+cont) = phi_num(cell_aux_v2(k+cont),1);
%                 localizacao_celulas_stencil(k+cont,:) = [cells(cell_aux_v2(k+cont),1) cells(cell_aux_v2(k+cont),2)];
%                 end 
if neuman && face_bound(stencil_faces(face_v2,cont),4) 
localizacao_celulas_stencil_neuman(cont,:) = [faces(stencil_faces(face_v2,cont),1) faces(stencil_faces(face_v2,cont),2) stencil_faces(face_v2,cont)];
end

  end 

%  factor=0;
%     if face_w_gauss && face_bound(face_v2,1)==1
%  solucao_numerica(guarda_face)=[];
%  solucao_analitica(guarda_face)=[];
%  factor=1;
%   end
 

tamanho21=size(localizacao_celulas_stencil);
plot(localizacao_celulas_stencil(:,1), localizacao_celulas_stencil(:,2),'*','color','r')




hold on
if neuman && ~isempty(localizacao_celulas_stencil_neuman)
plot(localizacao_celulas_stencil_neuman(:,1), localizacao_celulas_stencil_neuman(:,2),'*','color','g')
end
plot(xc, yc,'*','color','b')
 hold off    
axis([0 1 0 1])
legend(['FACE=' num2str(face_v2)]);
tamanho21(1)
[xc yc]
pause(0.1)
if tamanho21(1)==80 && xc>0.05
    pause(0.1)
end
 end  
 
 


end
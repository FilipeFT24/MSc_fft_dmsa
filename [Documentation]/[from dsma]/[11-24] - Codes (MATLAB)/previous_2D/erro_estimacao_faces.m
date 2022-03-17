close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       Filipe Diogo                                                %
%                                         07/03/2018                                               %
%                                                                                                   %
%                                                                                                   %
% Função que reconstroi os polnomios obtidos para cada face                                          %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



global L Lref cell_side vert_side cell_num face_num vert_num w;
global verts cells cell_verts cell_faces cell_vol faces face_area face_vert cell_norm;
global face_bound cell_bound face_cells vert_cells stencil_size stencil_cells stencil_faces T;
global phi lap_phi  phi_faces flux_phi_faces phi_num order solution fig face_w_gauss neuman T_border G constrained_source extended_stencil;

contador_kk=0;
mostra=[];

for face_v2 = 1:face_num
    
 erro_total(face_v2) = 0;
  erro_totalgradx(face_v2) =0;
   erro_totalgrady(face_v2) =0; 
   erro_total_c(face_v2) = 0;
  erro_totalgradx_c(face_v2) =0;
   erro_totalgrady_c(face_v2) =0; 
   
   num_erro_total(face_v2) = 0;
  num_erro_totalgradx(face_v2) =0;
   num_erro_totalgrady(face_v2) =0; 
   num_erro_total_c(face_v2) = 0;
  num_erro_totalgradx_c(face_v2) =0;
   num_erro_totalgrady_c(face_v2) =0; 

clear cell_aux_v2 solucao_analitica solucao_numerica localizacao_celulas_stencil C_v2 polinomio phi_estimado peso volume_celula volume_total
    
xf=faces(face_v2,1);
yf=faces(face_v2,2);

guarda_face=0;

 for k=1:stencil_size(face_v2,2)
                %
                % Celulas do Stencil %
                cell_aux_v2(k)=stencil_cells(face_v2,k);  
                solucao_analitica(k) = phi(1,cell_aux_v2(k));
                solucao_numerica(k) = phi_num(cell_aux_v2(k),1);
                localizacao_celulas_stencil(k,:) = [cells(cell_aux_v2(k),1) cells(cell_aux_v2(k),2)];  
 end    
 
  for cont=1:stencil_size(face_v2,3)
                %
                % faces do Stencil %
%                 if stencil_faces(face_v2,cont) ~= face_v2
                cell_aux_v2(k+cont)=stencil_faces(face_v2,cont);
                
%                 if neuman && face_bound(stencil_faces(face_v2,cont),4)
%                 solucao_numerica(k+cont) = SolutionDiffusion(solution,faces(face_v2,1),faces(face_v2,2),'xflux')*face_area(stencil_faces(face_v2,cont),1);
%                 solucao_analitica(k+cont) = SolutionDiffusion(solution,faces(face_v2,1),faces(face_v2,2),'xflux')*face_area(stencil_faces(face_v2,cont),1);
%                 else
                solucao_numerica(k+cont) = phi_faces(cell_aux_v2(k+cont));
                solucao_analitica(k+cont) = phi_faces(cell_aux_v2(k+cont));
%                 end
                localizacao_celulas_stencil(k+cont,:) = [faces(cell_aux_v2(k+cont),1) faces(cell_aux_v2(k+cont),2)];
if stencil_faces(face_v2,cont) == face_v2
 guarda_face = k+cont;
end
%                 else
%                 guarda_face = k+cont;
%                 cell_aux_v2(k+cont)=face_cells(face_v2,1); 
%                 solucao_analitica(k+cont) = phi(1,cell_aux_v2(k+cont));
%                 solucao_numerica(k+cont) = phi_num(cell_aux_v2(k+cont),1);
%                 localizacao_celulas_stencil(k+cont,:) = [cells(cell_aux_v2(k+cont),1) cells(cell_aux_v2(k+cont),2)];
%                 end                    
                
  end 

 factor=0;
 
%  if face_w_gauss && face_bound(face_v2,1)==1
%  solucao_numerica(guarda_face)=[];
%  solucao_analitica(guarda_face)=[];
%  factor=1;
%   end
 
% Constantes %
if face_w_gauss  && face_bound(face_v2,1)==1 


C_v2=T{face_v2}(:,:)*solucao_numerica';
C_v3=T{face_v2}(:,:)*solucao_analitica';

solucao_numerica(guarda_face)=[];
 solucao_analitica(guarda_face)=[];
 factor=1;
 
 
C_v2_c=T_border{face_v2}(:,:)*solucao_numerica'+constrained_source{face_v2};
C_v3_c=T_border{face_v2}(:,:)*solucao_analitica'+constrained_source{face_v2};
else
C_v2=T{face_v2}(:,:)*solucao_numerica';
C_v2_c=C_v2;
C_v3=T{face_v2}(:,:)*solucao_analitica';
C_v3_c=C_v3;
end


volume_total = 0; 
    
 
for qq=1:order/2
    
    xc=G{face_v2}(qq,1);
    yc=G{face_v2}(qq,2);     
     
if order ==2
phi_estimado(face_v2) = PolyReconstruction2ndOrder(C_v3,xc,yc,xf,yf,'poly',face_v2); 
phi_estimado_c(face_v2) = PolyReconstruction2ndOrder(C_v3_c,xc,yc,xf,yf,'poly',face_v2);

num_phi_estimado(face_v2) = PolyReconstruction2ndOrder(C_v2,xc,yc,xf,yf,'poly',face_v2); 
num_phi_estimado_c(face_v2) = PolyReconstruction2ndOrder(C_v2_c,xc,yc,xf,yf,'poly',face_v2);
elseif order == 4
phi_estimado(face_v2) = PolyReconstruction4thOrder(C_v3,xc,yc,xf,yf,'poly',face_v2);
phi_estimado_c(face_v2) = PolyReconstruction4thOrder(C_v3_c,xc,yc,xf,yf,'poly',face_v2);

num_phi_estimado(face_v2) = PolyReconstruction4thOrder(C_v2,xc,yc,xf,yf,'poly',face_v2);
num_phi_estimado_c(face_v2) = PolyReconstruction4thOrder(C_v2_c,xc,yc,xf,yf,'poly',face_v2);
elseif order ==6
phi_estimado(face_v2) = PolyReconstruction6thOrder(C_v3,xc,yc,xf,yf,'poly',face_v2);  
phi_estimado_c(face_v2) = PolyReconstruction6thOrder(C_v3_c,xc,yc,xf,yf,'poly',face_v2);

num_phi_estimado(face_v2) = PolyReconstruction6thOrder(C_v2,xc,yc,xf,yf,'poly',face_v2);  
num_phi_estimado_c(face_v2) = PolyReconstruction6thOrder(C_v2_c,xc,yc,xf,yf,'poly',face_v2);
elseif order ==8
phi_estimado(face_v2) = PolyReconstruction8thOrder(C_v3,xc,yc,xf,yf,'poly',face_v2);
phi_estimado_c(face_v2) = PolyReconstruction8thOrder(C_v3_c,xc,yc,xf,yf,'poly',face_v2);

num_phi_estimado(face_v2) = PolyReconstruction8thOrder(C_v2,xc,yc,xf,yf,'poly',face_v2);
num_phi_estimado_c(face_v2) = PolyReconstruction8thOrder(C_v2_c,xc,yc,xf,yf,'poly',face_v2);
else
   error('\n\nERRO: Reconstrução Não Implementada\n\n');
end


   
if order ==2
grad_phi_x(face_v2) = PolyReconstruction2ndOrder(C_v3,xc,yc,xf,yf,'xflux',face_v2); 
grad_phi_y(face_v2) = PolyReconstruction2ndOrder(C_v3,xc,yc,xf,yf,'yflux',face_v2);

grad_phi_x_c(face_v2) = PolyReconstruction2ndOrder(C_v3_c,xc,yc,xf,yf,'xflux',face_v2); 
grad_phi_y_c(face_v2) = PolyReconstruction2ndOrder(C_v3_c,xc,yc,xf,yf,'yflux',face_v2);

num_grad_phi_x(face_v2) = PolyReconstruction2ndOrder(C_v2,xc,yc,xf,yf,'xflux',face_v2); 
num_grad_phi_y(face_v2) = PolyReconstruction2ndOrder(C_v2,xc,yc,xf,yf,'yflux',face_v2);

num_grad_phi_x_c(face_v2) = PolyReconstruction2ndOrder(C_v2_c,xc,yc,xf,yf,'xflux',face_v2); 
num_grad_phi_y_c(face_v2) = PolyReconstruction2ndOrder(C_v2_c,xc,yc,xf,yf,'yflux',face_v2);
elseif order == 4
grad_phi_x(face_v2) = PolyReconstruction4thOrder(C_v3,xc,yc,xf,yf,'xflux',face_v2);  
grad_phi_y(face_v2) = PolyReconstruction4thOrder(C_v3,xc,yc,xf,yf,'yflux',face_v2);

grad_phi_x_c(face_v2) = PolyReconstruction4thOrder(C_v3_c,xc,yc,xf,yf,'xflux',face_v2);  
grad_phi_y_c(face_v2) = PolyReconstruction4thOrder(C_v3_c,xc,yc,xf,yf,'yflux',face_v2);


num_grad_phi_x(face_v2) = PolyReconstruction4thOrder(C_v2,xc,yc,xf,yf,'xflux',face_v2);  
num_grad_phi_y(face_v2) = PolyReconstruction4thOrder(C_v2,xc,yc,xf,yf,'yflux',face_v2);

num_grad_phi_x_c(face_v2) = PolyReconstruction4thOrder(C_v2_c,xc,yc,xf,yf,'xflux',face_v2);  
num_grad_phi_y_c(face_v2) = PolyReconstruction4thOrder(C_v2_c,xc,yc,xf,yf,'yflux',face_v2);

elseif order ==6
grad_phi_x(face_v2) = PolyReconstruction6thOrder(C_v3,xc,yc,xf,yf,'xflux',face_v2); 
grad_phi_y(face_v2) = PolyReconstruction6thOrder(C_v3,xc,yc,xf,yf,'yflux',face_v2);

grad_phi_x_c(face_v2) = PolyReconstruction6thOrder(C_v3_c,xc,yc,xf,yf,'xflux',face_v2); 
grad_phi_y_c(face_v2) = PolyReconstruction6thOrder(C_v3_c,xc,yc,xf,yf,'yflux',face_v2);

num_grad_phi_x(face_v2) = PolyReconstruction6thOrder(C_v2,xc,yc,xf,yf,'xflux',face_v2); 
num_grad_phi_y(face_v2) = PolyReconstruction6thOrder(C_v2,xc,yc,xf,yf,'yflux',face_v2);

num_grad_phi_x_c(face_v2) = PolyReconstruction6thOrder(C_v2_c,xc,yc,xf,yf,'xflux',face_v2); 
num_grad_phi_y_c(face_v2) = PolyReconstruction6thOrder(C_v2_c,xc,yc,xf,yf,'yflux',face_v2);

elseif order ==8
grad_phi_x(face_v2) = PolyReconstruction8thOrder(C_v3,xc,yc,xf,yf,'xflux',face_v2); 
grad_phi_y(face_v2) = PolyReconstruction8thOrder(C_v3,xc,yc,xf,yf,'yflux',face_v2);

grad_phi_x_c(face_v2) = PolyReconstruction8thOrder(C_v3_c,xc,yc,xf,yf,'xflux',face_v2); 
grad_phi_y_c(face_v2) = PolyReconstruction8thOrder(C_v3_c,xc,yc,xf,yf,'yflux',face_v2);

num_grad_phi_x(face_v2) = PolyReconstruction8thOrder(C_v2,xc,yc,xf,yf,'xflux',face_v2); 
num_grad_phi_y(face_v2) = PolyReconstruction8thOrder(C_v2,xc,yc,xf,yf,'yflux',face_v2);

num_grad_phi_x_c(face_v2) = PolyReconstruction8thOrder(C_v2_c,xc,yc,xf,yf,'xflux',face_v2); 
num_grad_phi_y_c(face_v2) = PolyReconstruction8thOrder(C_v2_c,xc,yc,xf,yf,'yflux',face_v2);

else
   error('\n\nERRO: Reconstrução Não Implementada\n\n');
end

analitico_phi(face_v2)=SolutionDiffusion(solution,xc,yc,'anal');
analitico_grad_x(face_v2)= SolutionDiffusion(solution,xc,yc,'xflux');
analitico_grad_y(face_v2)= SolutionDiffusion(solution,xc,yc,'yflux');

% if abs(phi_estimado(face_v2)-analitico_phi(face_v2)) >0.01
%     pause(1)
% end
  erro_total(face_v2) = erro_total(face_v2)+1/order*(abs(phi_estimado(face_v2)-analitico_phi(face_v2)));
  erro_total_c(face_v2) = erro_total_c(face_v2)+1/order*(abs(phi_estimado_c(face_v2)-analitico_phi(face_v2)));
  
  erro_totalgradx(face_v2) =erro_totalgradx(face_v2)+ 1/order*(abs(grad_phi_x(face_v2)-analitico_grad_x(face_v2)));
  erro_totalgradx_c(face_v2) =erro_totalgradx_c(face_v2)+ 1/order*(abs(grad_phi_x_c(face_v2)-analitico_grad_x(face_v2)));
  
  
   erro_totalgrady(face_v2) =erro_totalgrady(face_v2)+ 1/order*(abs(grad_phi_y(face_v2)-analitico_grad_y(face_v2))); 
   erro_totalgrady_c(face_v2) =erro_totalgrady_c(face_v2)+ 1/order*(abs(grad_phi_y_c(face_v2)-analitico_grad_y(face_v2))); 
   
   
    num_erro_total(face_v2) = num_erro_total(face_v2)+1/order*(abs(num_phi_estimado(face_v2)-analitico_phi(face_v2)));
  num_erro_total_c(face_v2) = num_erro_total_c(face_v2)+1/order*(abs(num_phi_estimado_c(face_v2)-analitico_phi(face_v2)));
  
  num_erro_totalgradx(face_v2) =num_erro_totalgradx(face_v2)+ 1/order*(abs(num_grad_phi_x(face_v2)-analitico_grad_x(face_v2)));
  num_erro_totalgradx_c(face_v2) =num_erro_totalgradx_c(face_v2)+ 1/order*(abs(num_grad_phi_x_c(face_v2)-analitico_grad_x(face_v2)));
  
  
   num_erro_totalgrady(face_v2) =num_erro_totalgrady(face_v2)+ 1/order*(abs(num_grad_phi_y(face_v2)-analitico_grad_y(face_v2))); 
   num_erro_totalgrady_c(face_v2) =num_erro_totalgrady_c(face_v2)+ 1/order*(abs(num_grad_phi_y_c(face_v2)-analitico_grad_y(face_v2))); 
   
end



if face_w_gauss  && face_bound(face_v2,1)==1 && yf==1
    contador_kk=contador_kk+1;
   erro_phi_only_border(contador_kk,:)= [xf yf erro_total(face_v2) erro_total_c(face_v2) num_erro_total(face_v2) num_erro_total_c(face_v2)];
   erro_gradxphi_only_border(contador_kk,:)= [xf yf erro_totalgradx(face_v2) erro_totalgradx_c(face_v2) num_erro_totalgradx(face_v2) num_erro_totalgradx_c(face_v2)];
   erro_gradyphi_only_border(contador_kk,:)= [xf yf erro_totalgrady(face_v2) erro_totalgrady_c(face_v2) num_erro_totalgrady(face_v2) num_erro_totalgrady_c(face_v2)];
end
 
 mostra(face_v2,:) = C_v3_c';
if face_bound(face_v2,1)==1 
[face_v2 face_bound(face_v2,1) ]  ;
mostra(face_v2,:);
pause(0.0001)
end
end

%   erro_total = erro_total';
%   erro_totalgradx = erro_totalgradx';
%    erro_totalgrady =  erro_totalgrady';
%    
%    media_phi=mean(erro_phi_only_border(:,3))
%    media_phi_c=mean(erro_phi_only_border(:,4))
%    num_media_phi=mean(erro_phi_only_border(:,5))
%    num_media_phi_c=mean(erro_phi_only_border(:,6))
%    
%    
%    media_gradx=mean(erro_gradxphi_only_border(:,3))
%    media_gradx_c=mean(erro_gradxphi_only_border(:,4))
%    num_media_gradx=mean(erro_gradxphi_only_border(:,5))
%    num_media_gradx_c=mean(erro_gradxphi_only_border(:,6))
%    
%    media_grady=mean(erro_gradyphi_only_border(:,3))
%    media_grady_C=mean(erro_gradyphi_only_border(:,4))
%    num_media_grady=mean(erro_gradyphi_only_border(:,5))
%    num_media_grady_C=mean(erro_gradyphi_only_border(:,6))
   

if face_w_gauss  && face_bound(face_v2,1)==1 
    
   max_phi=max(erro_phi_only_border(:,3))
   max_phi_c=max(erro_phi_only_border(:,4))
   num_max_phi=max(erro_phi_only_border(:,5))
   num_max_phi_c=max(erro_phi_only_border(:,6))
   
   max_gradx=max(erro_gradxphi_only_border(:,3))
    max_gradx_c=max(erro_gradxphi_only_border(:,4))
     num_max_gradx=max(erro_gradxphi_only_border(:,5))
    num_max_gradx_c=max(erro_gradxphi_only_border(:,6))
    
    
   max_grady=max(erro_gradyphi_only_border(:,3)) 
   max_grady_c=max(erro_gradyphi_only_border(:,4))
   num_max_grady=max(erro_gradyphi_only_border(:,5)) 
   num_max_grady_c=max(erro_gradyphi_only_border(:,6))

   
 
   media_phi_g = [mean(erro_phi_only_border(:,3)) mean(erro_phi_only_border(:,4)) mean(erro_phi_only_border(:,5)) mean(erro_phi_only_border(:,6))]
   media_gradx_g = [mean(erro_gradxphi_only_border(:,3)) mean(erro_gradxphi_only_border(:,4)) mean(erro_gradxphi_only_border(:,5)) mean(erro_gradxphi_only_border(:,6))]
   media_grady_g = [mean(erro_gradyphi_only_border(:,3)) mean(erro_gradyphi_only_border(:,4)) mean(erro_gradyphi_only_border(:,5)) mean(erro_gradyphi_only_border(:,6))]
   
   
end   
   
%    
%    soma_phi=sum(erro_total)
%    soma_gradx=sum(erro_totalgradx)
%    soma_grady=sum(erro_totalgrady)
% 
% k=0;
% for i=1:2:2*cell_side+1
%     for j=1:cell_side
%         k=k+1;
%         X(i,j)=faces(k,1);
%         Y(i,j)=faces(k,2);
%          modulo_errox(i,j)=erro_totalgradx(k);
%          modulo_erroy(i,j)=erro_totalgrady(k);
%         modulo_errophi(i,j)= erro_total(k);
%     end
% end
% for i=2:2:2*cell_side+1
%     for j=1:cell_side+1
%         k=k+1;
%         X(i,j)=faces(k,1);
%         Y(i,j)=faces(k,2);
%          modulo_errox(i,j)=erro_totalgradx(k);
%          modulo_erroy(i,j)=erro_totalgrady(k);
%         modulo_errophi(i,j)= erro_total(k);
%     end
% end
% 
% j=cell_side+1;
% for i=1:1:2*cell_side+1
%         if X(i,j)==0
%             X(i,j)=X(i,j-1);
%         end
%         if Y(i,j)==0
%             Y(i,j)=Y(i,j-1);
%         end
%         if modulo_errox(i,j)==0
%             modulo_errox(i,j)=modulo_errox(i,j-1);
%         end
%         if modulo_erroy(i,j)==0
%             modulo_erroy(i,j)=modulo_erroy(i,j-1);
%         end
%         if modulo_errophi(i,j)==0
%             modulo_errophi(i,j)=modulo_errophi(i,j-1);
%         end
% end
% % 
% % rr(cell_side+1,:)=[];
% % rr(1,:)=[];
% % rr(:,2*cell_side+1)=[];
% % rr(:,1)=[];
% % 
% % X(cell_side+1,:)=[];
% % X(1,:)=[];
% % X(:,2*cell_side+1)=[];
% % X(:,1)=[];
% % 
% % Y(cell_side+1,:)=[];
% % Y(1,:)=[];
% % Y(:,2*cell_side+1)=[];
% % Y(:,1)=[];
% % 
% fig(1)=figure;
% a=contourf(X,Y,modulo_errox);
% %legend('\phi_{analitico}');
% title('erro grad x face');
% xlabel('X');
% ylabel('Y');
% c=colorbar;
% 
% 
% fig(2)=figure;
% a=contourf(X,Y,modulo_erroy);
% %legend('\phi_{analitico}');
% title('erro grad y face');
% xlabel('X');
% ylabel('Y');
% c=colorbar;
% 
% fig(3)=figure;
% a=contourf(X,Y,modulo_errophi);
% %legend('\phi_{analitico}');
% title('erro phi face');
% xlabel('X');
% ylabel('Y');
% c=colorbar;

 pause(0.001)

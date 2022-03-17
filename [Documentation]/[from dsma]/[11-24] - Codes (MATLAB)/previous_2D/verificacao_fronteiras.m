
global L Lref cell_side vert_side cell_num face_num vert_num w;
global verts cells cell_verts cell_faces cell_vol faces face_area face_vert cell_norm;
global face_bound cell_bound face_cells vert_cells stencil_size stencil_cells stencil_faces T;
global phi lap_phi  phi_faces flux_phi_faces phi_num order fig T_border face_w_gauss constrained_source;

auxcc=0;
%  face_w_gauss=false;

for face_v2 = 1:face_num
    
    clear cell_aux_v2 solucao_analitica solucao_numerica localizacao_celulas_stencil C_v2 polinomio phi_estimado peso volume_celula volume_total
    
    if face_bound(face_v2,1)==1
        

         xf=faces(face_v2,1);
         yf=faces(face_v2,2);  


guarda_face=0;

% aux_yc(face_v2)= yc;
% aux_xc(face_v2)=xc;

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
               
                cell_aux_v2(k+cont)=stencil_faces(face_v2,cont);
                solucao_analitica(k+cont) = phi_faces(cell_aux_v2(k+cont));
                solucao_numerica(k+cont) = phi_faces(cell_aux_v2(k+cont));
                localizacao_celulas_stencil(k+cont,:) = [faces(cell_aux_v2(k+cont),1) faces(cell_aux_v2(k+cont),2)];
                
                if face_w_gauss
                   if stencil_faces(face_v2,cont) == face_v2
                     face_eli=k+cont;
                   end
                end
%                 else
%                 guarda_face = k+cont;
%                 cell_aux_v2(k+cont)=face_cells(face_v2,1); 
%                 solucao_analitica(k+cont) = phi(1,cell_aux_v2(k+cont));
%                 solucao_numerica(k+cont) = phi_num(cell_aux_v2(k+cont),1);
%                 localizacao_celulas_stencil(k+cont,:) = [cells(cell_aux_v2(k+cont),1) cells(cell_aux_v2(k+cont),2)];
%                 end                    
                
  end 

  if face_w_gauss
 solucao_numerica(face_eli)=[];
  end
% Constantes %

for qq=1:order/2
    
    xc=G{face_v2}(qq,1);
    yc=G{face_v2}(qq,2);
 
    
if face_w_gauss   
C_v2=T_border{face_v2}(:,:)*solucao_numerica'+constrained_source{face_v2};
else
C_v2=T{face_v2}(:,:)*solucao_numerica';
end

volume_total = 0; 
    
auxcc=auxcc+1;   
     
if order ==2
    if face_w_gauss   
phi_estimado_face(auxcc) = PolyReconstruction2ndOrder(C_v2,xc,yc,xf,yf,'poly');
grad_phi_x(auxcc) =  PolyReconstruction2ndOrder(C_v2,xc,yc,xf,yf,'xflux');
grad_phi_y(auxcc) =  PolyReconstruction2ndOrder(C_v2,xc,yc,xf,yf,'yflux');
  else
  phi_estimado_face(auxcc) = PolyReconstruction2ndOrder(C_v2,xc,yc,xf,yf,'poly');
  grad_phi_x(auxcc) =  PolyReconstruction2ndOrder(C_v2,xc,yc,xf,yf,'xflux');
  grad_phi_y(auxcc) =  PolyReconstruction2ndOrder(C_v2,xc,yc,xf,yf,'yflux');
  end
elseif order == 4
 
    if face_w_gauss   
phi_estimado_face(auxcc) = PolyReconstruction4thOrder(C_v2,xc,yc,xf,yf,'poly');
grad_phi_x(auxcc) =  PolyReconstruction4thOrder(C_v2,xc,yc,xf,yf,'xflux');
grad_phi_y(auxcc) =  PolyReconstruction4thOrder(C_v2,xc,yc,xf,yf,'yflux');
  else
  phi_estimado_face(auxcc) = PolyReconstruction4thOrder(C_v2,xc,yc,xf,yf,'poly');
  grad_phi_x(auxcc) =  PolyReconstruction4thOrder(C_v2,xc,yc,xf,yf,'xflux');
  grad_phi_y(auxcc) =  PolyReconstruction4thOrder(C_v2,xc,yc,xf,yf,'yflux');
  end
 
elseif order ==6
    
  if face_w_gauss   
phi_estimado_face(auxcc) = PolyReconstruction6thOrder(C_v2,xc,yc,xf,yf,'poly');
grad_phi_x(auxcc) =  PolyReconstruction6thOrder(C_v2,xc,yc,xf,yf,'xflux');
grad_phi_y(auxcc) =  PolyReconstruction6thOrder(C_v2,xc,yc,xf,yf,'yflux');
  else
  phi_estimado_face(auxcc) = PolyReconstruction6thOrder(C_v2,xc,yc,xf,yf,'poly');
  grad_phi_x(auxcc) =  PolyReconstruction6thOrder(C_v2,xc,yc,xf,yf,'xflux');
  grad_phi_y(auxcc) =  PolyReconstruction6thOrder(C_v2,xc,yc,xf,yf,'yflux');
  end


elseif order ==8
  if face_w_gauss   
phi_estimado_face(auxcc) = PolyReconstruction8thOrder(C_v2,xc,yc,xf,yf,'poly');
grad_phi_x(auxcc) =  PolyReconstruction8thOrder(C_v2,xc,yc,xf,yf,'xflux');
grad_phi_y(auxcc) =  PolyReconstruction8thOrder(C_v2,xc,yc,xf,yf,'yflux');
  else
  phi_estimado_face(auxcc) = PolyReconstruction8thOrder(C_v2,xc,yc,xf,yf,'poly');
  grad_phi_x(auxcc) =  PolyReconstruction8thOrder(C_v2,xc,yc,xf,yf,'xflux');
  grad_phi_y(auxcc) =  PolyReconstruction8thOrder(C_v2,xc,yc,xf,yf,'yflux');
  end
else
   error('\n\nERRO: Reconstrução Não Implementada\n\n');
end
 
analitico_phi(auxcc)=SolutionDiffusion('sin',xc,yc,'anal');
analitico_grad_x(auxcc)= SolutionDiffusion('sin',xc,yc,'xflux');
analitico_grad_y(auxcc)= SolutionDiffusion('sin',xc,yc,'yflux');

%  if abs(analitico_phi(auxcc)-phi_estimado_face(auxcc))>0.0000000000001
%      analitico_phi(auxcc)
%      phi_estimado_face(auxcc)
%      pause(0.1)
%  end

end

    end

end

 erro_medio_phi = mean(abs(phi_estimado_face-analitico_phi))
 erro_maximo_phi = max(abs(phi_estimado_face-analitico_phi))

 
 erro_medio_grad_x = mean(abs(grad_phi_x-analitico_grad_x))
 erro_maximo_grad_x = max(abs(grad_phi_x-analitico_grad_x))
 
  erro_medio_grad_y = mean(abs(grad_phi_y-analitico_grad_y))
 erro_maximo_grad_y = max(abs(grad_phi_y-analitico_grad_y))
 
figure()
plot(abs(phi_estimado_face-analitico_phi));
ylabel('erro')
xlabel('face_fronteira')
legend('Phi')

figure()
hold on
plot(abs(grad_phi_x-analitico_grad_x));
plot(abs(grad_phi_y-analitico_grad_y));
ylabel('erro')
xlabel('face_fronteira')
legend('Gradiente X', 'Gradiente Y')
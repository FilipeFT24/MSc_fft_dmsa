function [erro_medio, erro_maximo]=Polynomial_total_reconstruction()


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
global phi lap_phi  phi_faces flux_phi_faces phi_num order fig face_w_gauss T_border G constrained_source extended_stencil;


%  close(fig)

for face_v2 = 1:face_num
    
clear cell_aux_v2 solucao_analitica solucao_numerica localizacao_celulas_stencil C_v2 polinomio phi_estimado peso volume_celula volume_total
    
xc=faces(face_v2,1);
yc=faces(face_v2,2);

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
%                 if stencil_faces(face_v2,cont) ~= face_v2
                cell_aux_v2(k+cont)=stencil_faces(face_v2,cont);
                solucao_analitica(k+cont) = phi_faces(cell_aux_v2(k+cont));
                solucao_numerica(k+cont) = phi_faces(cell_aux_v2(k+cont));
                localizacao_celulas_stencil(k+cont,:) = [faces(cell_aux_v2(k+cont),1) faces(cell_aux_v2(k+cont),2)];
if stencil_faces(face_v2,cont) == face_v2
 guarda_face = k+cont;
 cell_aux_v2(k+cont)=face_cells(face_v2,1);
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
    if face_w_gauss && face_bound(face_v2,1)==1
 solucao_numerica(guarda_face)=[];
 solucao_analitica(guarda_face)=[];
 factor=1;
  end
 
% Constantes %
if face_w_gauss  && face_bound(face_v2,1)==1 
C_v2=T_border{face_v2}(:,:)*solucao_numerica'+constrained_source{face_v2};
else
C_v2=T{face_v2}(:,:)*solucao_numerica';
end


volume_total = 0; 
    
 for k=1:stencil_size(face_v2,2)
     
     
if order ==2
phi_estimado(k) = PolyReconstruction2ndOrder(C_v2,localizacao_celulas_stencil(k,1),localizacao_celulas_stencil(k,2),xc,yc,'poly',face_v2);    
elseif order == 4
phi_estimado(k) = PolyReconstruction4thOrder(C_v2,localizacao_celulas_stencil(k,1),localizacao_celulas_stencil(k,2),xc,yc,'poly',face_v2);   
elseif order ==6
phi_estimado(k) = PolyReconstruction6thOrder(C_v2,localizacao_celulas_stencil(k,1),localizacao_celulas_stencil(k,2),xc,yc,'poly',face_v2);    
elseif order ==8
phi_estimado(k) = PolyReconstruction8thOrder(C_v2,localizacao_celulas_stencil(k,1),localizacao_celulas_stencil(k,2),xc,yc,'poly',face_v2);
else
   error('\n\nERRO: Reconstrução Não Implementada\n\n');
end


dx=(cells(cell_aux_v2(k),1)-faces(face_v2,1));
dy=(cells(cell_aux_v2(k),2)-faces(face_v2,2));
peso(k)=(1/(sqrt(dx^2+dy^2))^order); 
volume_celula(k) = cell_vol(cell_aux_v2(k));
volume_total = volume_total+cell_vol(cell_aux_v2(k));

 end  
 
 for k=stencil_size(face_v2,2)+1:stencil_size(face_v2,2)+stencil_size(face_v2,3)-1*factor
     
     
if order ==2
phi_estimado(k) = PolyReconstruction2ndOrder(C_v2,localizacao_celulas_stencil(k,1),localizacao_celulas_stencil(k,2),xc,yc,'poly',face_v2);    
elseif order == 4
phi_estimado(k) = PolyReconstruction4thOrder(C_v2,localizacao_celulas_stencil(k,1),localizacao_celulas_stencil(k,2),xc,yc,'poly',face_v2);   
elseif order ==6
phi_estimado(k) = PolyReconstruction6thOrder(C_v2,localizacao_celulas_stencil(k,1),localizacao_celulas_stencil(k,2),xc,yc,'poly',face_v2);    
elseif order ==8
phi_estimado(k) = PolyReconstruction8thOrder(C_v2,localizacao_celulas_stencil(k,1),localizacao_celulas_stencil(k,2),xc,yc,'poly',face_v2);
else
   error('\n\nERRO: Reconstrução Não Implementada\n\n');
end


dx=(faces(cell_aux_v2(k),1)-faces(face_v2,1));
dy=(faces(cell_aux_v2(k),2)-faces(face_v2,2));
if k== guarda_face
dx=(cells(cell_aux_v2(k),1)-faces(face_v2,1));
dy=(cells(cell_aux_v2(k),2)-faces(face_v2,2));
%     dx=  0.0000001;
%     dy = 0.0000001;
peso(k)=(1/(sqrt(dx^2+dy^2))^order); 
else
peso(k)=(1/(sqrt(dx^2+dy^2))^order); 
end

if peso(k)==Inf
   pause(0.1) 
end


 end  
     
erro_total(face_v2) = mean((abs(phi_estimado-solucao_analitica)'));
erro_total2(face_v2) = mean((peso'.*abs(phi_estimado-solucao_numerica)')/sum(peso));
SSE(face_v2) =  sum(((peso'.*abs(phi_estimado-solucao_numerica)')/sum(peso)).^2);
SSTO(face_v2) =  sum(((peso'.*abs(mean(solucao_numerica)^2-solucao_numerica.^2)')/sum(peso)));
% SSTO(face_v2) =  sum((abs(mean(solucao_numerica)-solucao_numerica)').^2);

r_squared(face_v2)=1-SSE(face_v2)/SSTO(face_v2);
r_squared_minusone(face_v2)=SSE(face_v2)/SSTO(face_v2);
pontos_por_face(face_v2)=stencil_size(face_v2,2)+stencil_size(face_v2,3)-1*factor;
 
% figure()
% plot(1./sqrt(peso),abs(phi_estimado-solucao_analitica),'*', 'color', 'b')
% figure()
% plot(1./sqrt(peso),(peso'.*abs(phi_estimado-solucao_analitica)')/sum(peso),'*', 'color', 'r') 
end

% figure()
% plot3(aux_xc,aux_yc,erro_total,'*')

erro_medio = mean(erro_total2);
erro_maximo = max(erro_total2);


k=0;
for j=1:2:2*cell_side+1
    for i=1:cell_side
        k=k+1;
        X(i,j)=faces(k,1);
        Y(i,j)=faces(k,2);
        rr(i,j)=r_squared(k);
        rr2(i,j)=r_squared_minusone(k);
        modulo_erro(i,j)= erro_total(k);
        erro2(i,j)= erro_total2(k);
        num_points(i,j)= pontos_por_face(k);
    end
end
for j=2:2:2*cell_side+1
    for i=1:cell_side+1
        k=k+1;
        X(i,j)=faces(k,1);
        Y(i,j)=faces(k,2);
        rr(i,j)=r_squared(k);
        rr2(i,j)=r_squared_minusone(k);
        modulo_erro(i,j)= erro_total(k);
        erro2(i,j)= erro_total2(k);
        num_points(i,j)= pontos_por_face(k);
    end
end

i=cell_side+1;
for j=1:1:2*cell_side+1
        if X(i,j)==0
            X(i,j)=X(i-1,j);
        end
        if Y(i,j)==0
            Y(i,j)=Y(i-1,j);
        end
        if rr(i,j)==0
            rr(i,j)=rr(i-1,j);
        end
        if rr2(i,j)==0
            rr2(i,j)=rr2(i-1,j);
        end
        if modulo_erro(i,j)==0
            modulo_erro(i,j)=modulo_erro(i-1,j);
        end
        if erro2(i,j)==0
            erro2(i,j)=erro2(i-1,j);
        end
          if num_points(i,j)==0
            num_points(i,j)=num_points(i-1,j);
        end
end
% 
rr(cell_side+1,:)=[];
rr(1,:)=[];
rr(:,2*cell_side+1)=[];
rr(:,1)=[];

X(cell_side+1,:)=[];
X(1,:)=[];
X(:,2*cell_side+1)=[];
X(:,1)=[];

Y(cell_side+1,:)=[];
Y(1,:)=[];
Y(:,2*cell_side+1)=[];
Y(:,1)=[];

rr2(cell_side+1,:)=[];
rr2(1,:)=[];
rr2(:,2*cell_side+1)=[];
rr2(:,1)=[];

modulo_erro(cell_side+1,:)=[];
modulo_erro(1,:)=[];
modulo_erro(:,2*cell_side+1)=[];
modulo_erro(:,1)=[];

erro2(cell_side+1,:)=[];
erro2(1,:)=[];
erro2(:,2*cell_side+1)=[];
erro2(:,1)=[];


num_points(cell_side+1,:)=[];
num_points(1,:)=[];
num_points(:,2*cell_side+1)=[];
num_points(:,1)=[];

% fig(1)=figure;
% a=contourf(X,Y,rr);
% %legend('\phi_{analitico}');
% title('R^2');
% xlabel('X');
% ylabel('Y');
% c=colorbar;
% 
% fig(4)=figure;
% a=contourf(X,Y,rr2);
% %legend('\phi_{analitico}');
% title('1-R^2');
% xlabel('X');
% ylabel('Y');
% c=colorbar;
% 
% fig(2)=figure;
% a=contourf(X,Y,modulo_erro);
% %legend('\phi_{analitico}');
% title('erro (Analitico vs Polinomio)');
% xlabel('X');
% ylabel('Y');
% c=colorbar;
% 
% fig(3)=figure;
% a=contourf(X,Y,erro2);
% %legend('\phi_{analitico}');
% title('erro (numerico vs Polinomio)');
% xlabel('X');
% ylabel('Y');
% c=colorbar;
% 
% 
% % fig(3)=figure;
% figure()
% a=contourf(X,Y,num_points);
% %legend('\phi_{analitico}');
% title('Nº Pontos do Stencil');
% xlabel('X');
% ylabel('Y');
% c=colorbar;


% pause(0.001)
end
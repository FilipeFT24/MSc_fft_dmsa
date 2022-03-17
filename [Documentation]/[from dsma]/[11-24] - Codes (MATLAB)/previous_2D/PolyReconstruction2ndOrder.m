function f=PolyReconstruction2ndOrder(C,x,y,xc,yc,type,face)

global L Lref cell_side vert_side cell_num face_num vert_num w;
global verts cells cell_verts cell_faces cell_vol faces face_area face_vert cell_norm;
global face_bound cell_bound face_cells vert_cells;
global phi lap_phi  phi_faces flux_phi_faces face_w_gauss extended_stencil;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                  14 de Outubro de 2016                                            %
%                                  14 de Outubro de 2016                                            %
%                                                                                                   %
% Fun��o calcula a contirbui��o de cada celula para a reconstru��o do polinomio na face que se esta %
% a trabalhar, esta fun��o tem como objectivo facilitar a constru��o da matriz A.                   %
% Aten��o que esta fun��o n�o determina o polinomio reconstruido na face, serve apenas para quando  %
% se esta a determinar as entradas na matriz A.                                                     %
% Caso queira saber o polinomio da face tenho de usar outra fun��o                                  %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    % Normal Exterior da Face %
    vert=face_vert(face,:);
    nx=verts(vert(1),1)-verts(vert(2),1);
    ny=verts(vert(1),2)-verts(vert(2),2);



if strcmp(type,'poly')==1
    if extended_stencil && nx==0
    f=C(1)+C(2)*(x-xc)+C(3)*(y-yc) ...
        +C(4)*(y-yc)^2+C(5)*(y-yc)*(x-xc);
    elseif extended_stencil && ny==0
    f=C(1)+C(2)*(x-xc)+C(3)*(y-yc) ...
    +C(4)*(x-xc)^2+C(5)*(y-yc)*(x-xc);
    else
    f=C(1)+C(2)*(x-xc)+C(3)*(y-yc);
%     C(1)
%     C(2)
%     C(3)
%     x-xc
%     y-yc
%     pause(0.1)
    end
elseif strcmp(type,'xflux')==1
     if extended_stencil && nx==0
         f=C(2) ...
             +C(5)*(y-yc);
    elseif extended_stencil && ny==0
        f=C(2) ...
            +2*C(4)*(x-xc)+C(5)*(y-yc);
     else
    f=C(2);
     end
elseif strcmp(type,'yflux')==1
     if extended_stencil && nx==0
         f=C(3) ...
             +2*C(4)*(y-yc)+C(5)*(x-xc);
    elseif extended_stencil && ny==0
        f=C(3) ...
            +C(5)*(x-xc);
     else
    f=C(3);
     end
else
    error('\n\n\tERRO: Fluxo Desconhecido\n\n');
end
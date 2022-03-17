function [G]=GaussFace(order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                  11 de Outubro de 2016                                            %
%                                  13 de Outubro de 2016                                            %
%                                                                                                   %
%                                                                                                   %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global L Lref cell_side vert_side cell_num face_num vert_num;
global verts cells cell_verts cell_faces cell_vol faces face_area face_vert;
global face_bound cell_bound face_cells vert_cells;
%
type='1D';
%
for i=1:face_num
    %
    % Face Selecionada %
    face=i;
    %
    % Vertices da Face %
    vert1=face_vert(face,1);
    vert2=face_vert(face,2);
    %
    % Coordenadas do Vertice 1 %
    x1=verts(vert1,1);
    y1=verts(vert1,2);
    %
    % Coordenadas do Vertice 2 %
    x2=verts(vert2,1);
    y2=verts(vert2,2);
    %
    % Coordenadas do Centroide da Face %
    x3=faces(face,1);
    y3=faces(face,2);
    %
    % Calcula as Coordenadas dos Pontos de Gauss %
    if order==2
        G{face}(1,1)=x3;
        G{face}(1,2)=y3;
        G{face}(1,3)=1;
    else
        [G{face},n]=gausspoints(x1,y1,x2,y2,x3,y3,type,order);
    end
end
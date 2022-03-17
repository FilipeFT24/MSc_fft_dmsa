function [pontos,n]=gausspoints(x1,y1,x2,y2,x3,y3,type,order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   28 de Junho de 2016                                             %
%                                                                                                   %
% Função que determina as soluções analiticas                                                       %
%                                                                                                   %
% x1        - Coordenadas do Vertice 1                                                              %
% x2        - Coordenadas do Vertice 2                                                              %
% x3        - Coordenadas do Centroide da Face ou da Célula                                         %
% type      - Tipo de Integração que se pretende fazer, 1D ou 2D, sendo uma string                  %
% order     - Ordem da Integração que se está a realizar                                            %
%                                                                                                   %
% pontos    - Coordenadas dos Pontos de Gauss,(x,y,peso), as coordenadas são reais, ou seja já se   %
%             tem em conta a posição do centroide da face ou da célula que se está a fazer os       %
%             cálculos                                                                              %
% n         - Numero de Pontos de Gauss;                                                            %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global L Lref cell_side vert_side cell_num face_num vert_num;
global verts cells cell_verts cell_faces cell_vol faces face_area face_vert;
global face_bound cell_bound face_cells vert_cells;
%
%% 1D %%
%
if type=='1D'
    dx=abs(x2-x1)/2;
    dy=abs(y2-y1)/2;
    %
    % 4ª Ordem %
    if order==4
        %
        % Ponto 1 %
        pontos(1,1)=x3+dx*sqrt(1/3);
        pontos(1,2)=y3+dy*sqrt(1/3);
        pontos(1,3)=1/2;
        %
        % Ponto 2 %
        pontos(2,1)=x3-dx*sqrt(1/3);
        pontos(2,2)=y3-dy*sqrt(1/3);
        pontos(2,3)=1/2;
        %
        n=2;
        %
        % 6ª Ordem % 
    elseif order==6
        %
        % Ponto 1 %
        pontos(1,1)=x3+dx*sqrt(3/5);
        pontos(1,2)=y3+dy*sqrt(3/5);
        pontos(1,3)=5/18;
        % 
        % Ponto 2 %
        pontos(2,1)=x3;
        pontos(2,2)=y3;
        pontos(2,3)=4/9;
        % 
        % Ponto 3 %
        pontos(3,1)=x3-dx*sqrt(3/5);
        pontos(3,2)=y3-dy*sqrt(3/5);
        pontos(3,3)=5/18;
        %
        n=3;
        %
        % 8ª Ordem %
    elseif order==8
        %
        % Ponto 1 %
        pontos(1,1)=x3+dx*sqrt(525+70*sqrt(30))/35;
        pontos(1,2)=y3+dy*sqrt(525+70*sqrt(30))/35;
        pontos(1,3)=(18-sqrt(30))/72;
        %
        % Ponto 2 %
        pontos(2,1)=x3+dx*sqrt(525-70*sqrt(30))/35;
        pontos(2,2)=y3+dy*sqrt(525-70*sqrt(30))/35;
        pontos(2,3)=(18+sqrt(30))/72;
        % 
        % Ponto 3 %
        pontos(3,1)=x3-dx*sqrt(525-70*sqrt(30))/35;
        pontos(3,2)=y3-dy*sqrt(525-70*sqrt(30))/35;
        pontos(3,3)=(18+sqrt(30))/72;
        % 
        % Ponto 4 %
        pontos(4,1)=x3-dx*sqrt(525+70*sqrt(30))/35;
        pontos(4,2)=y3-dy*sqrt(525+70*sqrt(30))/35;
        pontos(4,3)=(18-sqrt(30))/72;
        %
        n=4;
    else
        error('\n\n\tERRO: Método Desconhecido\n\n');
    end
    %
    %% 2D %%
    %
elseif type=='2D'
    %
    % 4ª Ordem %
    if order==4
        %
        % Ponto 1 %
        pontos(1,1)=(x1+x2+x3)/3;
        pontos(1,2)=(y1+y2+y3)/3;
        pontos(1,3)=-9/16;
        %
        % Ponto 2 %
        pontos(2,1)=x1*3/5+(x2+x3)/5;
        pontos(2,2)=y1*3/5+(y2+y3)/5;
        pontos(2,3)=25/48;
        % 
        % Ponto 3 %
        pontos(3,1)=x2*3/5+(x1+x3)/5;
        pontos(3,2)=y2*3/5+(y1+y3)/5;
        pontos(3,3)=25/48;
        % 
        % Ponto 4 %
        pontos(4,1)=x3*3/5+(x1+x2)/5;
        pontos(4,2)=y3*3/5+(y1+y2)/5;
        pontos(4,3)=25/48;
        %
        n=4;
        %
        % 6ª Ordem %
    elseif order==6
        %
        % Constantes %
        a1=15003009/18814273;
        b1=1905632/18814273;
        a2=952816/15955825;
        b2=15003009/31911650;
        %
        % Ponto 1 %
        pontos(1,1)=(x1+x2+x3)/3;
        pontos(1,2)=(y1+y2+y3)/3;
        pontos(1,3)=9/40;
        %
        % Ponto 2 %
        pontos(2,1)=x1*a1+(x2+x3)*b1;
        pontos(2,2)=y1*a1+(y2+y3)*b1;
        pontos(2,3)=4057154/32215185;
        % 
        % Ponto 3 %
        pontos(3,1)=x2*a1+(x1+x3)*b1;
        pontos(3,2)=y2*a1+(y1+y3)*b1;
        pontos(3,3)=4057154/32215185;
        % 
        % Ponto 4 %
        pontos(4,1)=x3*a1+(x1+x2)*b1;
        pontos(4,2)=y3*a1+(y1+y2)*b1;
        pontos(4,3)=4057154/32215185;
        % 
        % Ponto 5 %
        pontos(5,1)=x1*a2+(x2+x3)*b2;
        pontos(5,2)=y1*a2+(y2+y3)*b2;
        pontos(5,3)=8846073/66816191;
        % 
        % Ponto 6 %
        pontos(6,1)=x2*a2+(x1+x3)*b2;
        pontos(6,2)=y2*a2+(y1+y3)*b2;
        pontos(6,3)=8846073/66816191;
        % 
        % Ponto 7 %
        pontos(7,1)=x3*a2+(x1+x2)*b2;
        pontos(7,2)=y3*a2+(y1+y2)*b2;
        pontos(7,3)=8846073/66816191;
        %
        n=7;
        %
        % 8ª Ordem %
    elseif order==8
        %
        % Contates %
        a1=31723247/66185506;
        b1=10026259/38511290;
        a2=34357532/39503231;
        b2=2273498/34907023;
        a3=990949/20352076;
        b3=7686421/24567813;
        c3=6574983/10298446;
        % 
        % Ponto 1 %
        pontos(1,1)=(x1+x2+x3)/3;
        pontos(1,2)=(y1+y2+y3)/3;
        pontos(1,3)=-3294447/22026115;
        % 
        % Ponto 2 %
        pontos(2,1)=x1*a1+(x2+x3)*b1;
        pontos(2,2)=y1*a1+(y2+y3)*b1;
        pontos(2,3)=21277533/121159934;
        % 
        % Ponto 3 %
        pontos(3,1)=x2*a1+(x1+x3)*b1;
        pontos(3,2)=y2*a1+(y1+y3)*b1;
        pontos(3,3)=21277533/121159934;
        % 
        % Ponto 4 %
        pontos(4,1)=x3*a1+(x1+x2)*b1;
        pontos(4,2)=y3*a1+(y1+y2)*b1;
        pontos(4,3)=21277533/121159934;
        % 
        % Ponto 5 %
        pontos(5,1)=x1*a2+(x2+x3)*b2;
        pontos(5,2)=y1*a2+(y2+y3)*b2;
        pontos(5,3)=2259506/42354697;
        % 
        % Ponto 6 %
        pontos(6,1)=x2*a2+(x1+x3)*b2;
        pontos(6,2)=y2*a2+(y1+y3)*b2;
        pontos(6,3)=2259506/42354697;
        % 
        % Ponto 7 %
        pontos(7,1)=x3*a2+(x1+x2)*b2;
        pontos(7,2)=y3*a2+(y1+y2)*b2;
        pontos(7,3)=2259506/42354697;
        % 
        % Ponto 8 %
        pontos(8,1)=x1*a3+x2*b3+x3*c3;
        pontos(8,2)=y1*a3+y2*b3+y3*c3;
        pontos(8,3)=8113102/105209523;
        % 
        % Ponto 9 %
        pontos(9,1)=x1*a3+x3*b3+x2*c3;
        pontos(9,2)=y1*a3+y3*b3+y2*c3;
        pontos(9,3)=8113102/105209523;
        % 
        % Ponto 10 %
        pontos(10,1)=x2*a3+x1*b3+x3*c3;
        pontos(10,2)=y2*a3+y1*b3+y3*c3;
        pontos(10,3)=8113102/105209523;
        % 
        % Ponto 11 %
        pontos(11,1)=x2*a3+x3*b3+x1*c3;
        pontos(11,2)=y2*a3+y3*b3+y1*c3;
        pontos(11,3)=8113102/105209523;
        % 
        % Ponto 12 %
        pontos(12,1)=x3*a3+x1*b3+x2*c3;
        pontos(12,2)=y3*a3+y1*b3+y2*c3;
        pontos(12,3)=8113102/105209523;
        % 
        % Ponto 12 %
        pontos(13,1)=x3*a3+x2*b3+x1*c3;
        pontos(13,2)=y3*a3+y2*b3+y1*c3;
        pontos(13,3)=8113102/105209523;
        %
        n=13;
    end
end
%
end

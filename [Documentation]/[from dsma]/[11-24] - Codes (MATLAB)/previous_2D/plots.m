function a=plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   04 de Julho de 2016                                             %
%                                                                                                   %
% Função que cria os graficos do problema                                                           %
%                                                                                                   %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global L Lref cell_side vert_side cell_num face_num vert_num;
global verts cells cell_verts cell_faces cell_vol faces face_area face_vert;
global face_bound cell_bound face_cells vert_cells;
global phi lap_phi  phi_faces flux_phi_faces;
global phi_num lap_phi_num;
global norma1_phi erro_phi_max erro_phi norma1_lap erro_lap_max erro_lap;
%
%% Variaveis Auxiliares %%
% Malha %
k=0;
for i=1:cell_side
    for j=1:cell_side
        k=k+1;
        X(i,j)=cells(k,1);
        Y(i,j)=cells(k,2);
    end
end
% Valores Analiticos %
k=0;
for i=1:cell_side
    for j=1:cell_side
        k=k+1;
        Z1(i,j)=phi(k);
        Z2(i,j)=lap_phi(k);
        erro_maximo(i,j)=erro_phi(k);
        result_num(i,j)=phi_num(k);
    end
end
%% Plots %%
%
% Phi Analitico %
figure (1)
a=contourf(X,Y,Z1);
%legend('\phi_{analitico}');
title('\phi_{analitico}');
xlabel('X');
ylabel('Y');
c=colorbar;
c.Label.String = '\phi_{analitico}';
% Phi Analitico %
figure (6)
a=contourf(X,Y,result_num);
%legend('\phi_{analitico}');
title('\phi_{numerico}');
xlabel('X');
ylabel('Y');
c=colorbar;

% Laplaciano Analitico %
% figure (2)
% contourf(X,Y,Z2)
% %legend('\nabla^2\phi_{analitico}');
% title('\nabla^2\phi_{analitico}');
% xlabel('X');
% ylabel('Y');
% colorbar;

figure (3)
contourf(X,Y,erro_maximo)
% ERRO phi max
title('\phi_{analitico} Erro');
xlabel('X');
ylabel('Y');
colorbar;
figure (4)
surf(X,Y,Z1)
% ERRO phi max
title('\phi_{analitico}');
xlabel('X');
ylabel('Y');
colorbar;

figure (5)
surf(X,Y,result_num)
% ERRO phi max
title('\phi_{numerico}');
xlabel('X');
ylabel('Y');
colorbar;


figure (6)
surf(X,Y,result_num)
hold on
surf(X,Y,Z1)
% ERRO phi max
title('\phi_{numerico} vs analitico');
xlabel('X');
ylabel('Y');
colorbar;

%

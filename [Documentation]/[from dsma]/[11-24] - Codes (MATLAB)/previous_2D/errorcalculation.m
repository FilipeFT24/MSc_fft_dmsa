function [norma1_phi,norma1_lap_phi,erro_max_phi,erro_max_lap_phi,erro_phi,erro_lap_phi,X]=errorcalculation(explicito)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   22 de Setembro de 2016                                          %
%                                                                                                   %
% Cálculo dos Erros Numéricos                                                                       %
%                                                                                                   %
% norma1_phi - Norma 1 do Erro da Propriedade                                                       %
% norma1_lap - Norma 1 do Erro do Laplaciano                                                        %
% erro_phi - Distribuição do Erro da Propriedade                                                    %
% erro_lap - Distribuição do Erro do Laplaciano                                                     %
% erro_phi_max - Valor do Erro Máximo da Propriedade                                                %
% erro_lap_max - Valor do Erro Máximo do Laplaciano                                                 &
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Declaração das Variaveis Globais %%
%
global L Lref cell_side vert_side cell_num face_num vert_num;
global verts cells cell_verts cell_faces cell_vol faces face_area face_vert;
global face_bound cell_bound face_cells vert_cells;
global phi lap_phi  phi_faces flux_phi_faces;
global phi_num lap_phi_num;
%
%% Cálculo do Erro %%
%
erro_max_phi=0;
erro_max_lap_phi=0;
erro_phi=0;
erro_lap_phi=0;
erro_phi_total=0;
erro_lap_phi_total=0;
%
if explicito
    % Nada %
else
    for i=1:cell_num
        erro_lap_phi(i)=abs(lap_phi(i)-lap_phi_num(i));
        erro_phi(i)=abs(phi_num(i)-phi(i));
        erro_lap_phi_total=erro_lap_phi_total+erro_lap_phi(i)*cell_vol(i);
        erro_phi_total=erro_phi_total+erro_phi(i)*cell_vol(i);
        %
        if erro_max_lap_phi<erro_lap_phi(i)
            erro_max_lap_phi=erro_lap_phi(i);
            X(1,1)=cells(i,1);
            X(1,2)=cells(i,2);
        end
        if erro_max_phi<erro_phi(i)
            erro_max_phi=erro_phi(i);
            X(2,1)=cells(i,1);
            X(2,2)=cells(i,2);
        end
    end
end
%
norma1_lap_phi=erro_lap_phi_total/L^2;
norma1_phi=erro_phi_total/L^2;
%
end
function [norma1_phi,erro_phi_max,erro_phi,norma1_lap,erro_lap_max,erro_lap, aux4, aux2]=errorcalculation(explicito)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   07 de Junho de 2016                                             %
%                                                                                                   %
% Cálculo dos Erros Numéricos                                                                       %
%                                                                                                   %
% norma1_phi - Norma 1 do Erro da Propriedade                                                       %
% norma1_lap - Norma 1 do Erro do Laplaciano                                                        %
% erro_phi - Distribuição do Erro da Propriedade                                                    %
% erro_lap - Distribuição do Erro do Laplaciano                                                     %
% erro_phi_max - Valor do Erro Máximo da Propriedade                                                %
% erro_lap_max - Valor do Erro Máximo do Laplaciano                                                 &
% aux2 - Norma 1 do Erro sem a fronteira                                                            %
% aux4 - Erro Máximo sem a fronteira                                                                %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Declaração das Variaveis Globais %%
%
global L Lref cell_num face_num;
global x faces cell_vol;
global phi lap_phi phi_Weast phi_East;
global phi_num lap_phi_num;
%
%% Cálculo do Erro %%
%
erro_phi_max=0;
erro_phi_total=0;
erro_lap_max=0;
erro_lap_total=0;
aux1=0;
aux2=0;
aux3=0;
%
if explicito
    for i=1:cell_num
        erro_lap(i)=abs(lap_phi(i)-lap_phi_num(i));
        erro_phi(i)=0;
        erro_lap_total=erro_lap_total+erro_lap(i)*cell_vol(i);
        if erro_lap_max<erro_lap(i)
            erro_lap_max=erro_lap(i);
        end
    end
    for i=4:cell_num-3
        aux1=aux1+erro_lap(i)*cell_vol(i);
        if aux2<erro_lap(i)
            aux2=erro_lap(i);
        end
        aux3=aux3+cell_vol(i);
    end
else
    for i=1:cell_num
        erro_lap(i)=abs(lap_phi(i)-lap_phi_num(i));
        erro_lap_total=erro_lap_total+erro_lap(i)*cell_vol(i);
        erro_phi(i)=abs(phi(i)-phi_num(i));
        erro_phi_total=erro_phi_total+erro_phi(i)*cell_vol(i);
        if erro_phi_max<erro_phi(i)
            erro_phi_max=erro_phi(i);
        end
        if erro_lap_max<erro_lap(i)
            erro_lap_max=erro_lap(i);
        end
    end
    for i=4:cell_num-3
        aux1=aux1+erro_phi(i)*cell_vol(i);
        if aux2<erro_phi(i)
            aux2=erro_phi(i);
        end
        aux3=aux3+cell_vol(i);
    end
    
end
%
norma1_lap=erro_lap_total/L;
norma1_phi=erro_phi_total/L;
%
aux4=aux1/aux3;
%
end
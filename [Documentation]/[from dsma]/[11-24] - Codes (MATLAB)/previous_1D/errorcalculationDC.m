function [norma1_phi,erro_phi_max,erro_phi]=errorcalculationDC
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
%
for i=1:cell_num
    
    erro_phi(i)=abs(phi(i)-phi_num(i));
    erro_phi_total=erro_phi_total+erro_phi(i)*cell_vol(i);
    if erro_phi_max<erro_phi(i)
        erro_phi_max=erro_phi(i);
    end
end
%
norma1_phi=erro_phi_total/L;
%
end
function [phi_num,lap_phi_num]=Scheme2ndOrder(explicito,dirichlet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   07 de Junho de 2016                                             %
%                                   27 de Junho de 2016                                             %
%                                                                                                   %
% Resolução do Problema de 2ª Ordem                                                                 %
%                                                                                                   %
% Condições de Fronteira de Dirichlet nas duas fronteiras;                                          %
% Condições de Fronteira de Neumann na fronteira esquerda e de Dirichlet na fronteira esquerda;     %
%                                                                                                   %
% Cálculo Explicito é a verificação dos fluxos, reconstroi-se o polinomio nas faces com base em     %
% diferenças finitas de 2ª ordem e depois a partir do método de volume finito calcula-se os fluxos  %
% nas faces da célula.                                                                              %
%                                                                                                   %
% Cálculo Implicito é a determinação da propriedade no centroide da célula com base nas             %
% reconstruçoes realizadas para as faces, sendo utilizado o método de volume finito, neste situação %
% é utilizado o algoritmo de Thomas para se obter o resultado final.                                %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Declaração das Variaveis Globais %%
%
global TempoInverterMatriz;
global L Lref cell_num face_num;
global x faces cell_vol;
global phi lap_phi phi_West phi_East flux_phi_West flux_phi_East;
%
%%
%
if explicito
    TempoInverterMatriz=cputime;
    for i=1:cell_num
        if i==1
            if dirichlet
                flux_phi_west=(phi(i)-phi_West)/x(i);
                flux_phi_east=(phi(i+1)-phi(i))/(x(i+1)-x(i));
            else
                flux_phi_west=flux_phi_West;
                flux_phi_east=(phi(i+1)-phi(i))/(x(i+1)-x(i));
            end
        elseif i==cell_num
            flux_phi_west=(phi(i)-phi(i-1))/(x(i)-x(i-1));
            flux_phi_east=(phi_East-phi(i-1))/(L-x(i));
        else 
            flux_phi_west=(phi(i)-phi(i-1))/(x(i)-x(i-1));
            flux_phi_east=(phi(i+1)-phi(i))/(x(i+1)-x(i));
        end
        lap_phi_num(i)=(flux_phi_east-flux_phi_west)/cell_vol(i);
        phi_num(i)=phi(i);
    end
    TempoInverterMatriz=cputime-TempoInverterMatriz;
else
    TempoInverterMatriz=cputime;
    
    
    % Algoritmo Thomas - Condensação %
    A=zeros(cell_num);
    for i=1:cell_num
        if i==1
            if dirichlet
                a(i)=0;
                d(i)=-(1/(x(i+1)-x(i))+1/x(i));
                c(i)=1/(x(i+1)-x(i));
                b(i)=lap_phi(i)*cell_vol(i)-phi_West/x(i);
            else
                a(i)=0;
                d(i)=-1/(x(i+1)-x(i));
                c(i)=1/(x(i+1)-x(i));
                b(i)=lap_phi(i)*cell_vol(i)+flux_phi_West;
            end
            %
            d1(i)=d(i);
            b1(i)=b(i);
        elseif i==cell_num
            a(i)=1/(x(i)-x(i-1));
            d(i)=-(1/(L-x(i))+1/(x(i)-x(i-1)));
            c(i)=0;
            b(i)=lap_phi(i)*cell_vol(i)-phi_East/(L-x(i));
            %
            d1(i)=d(i)-a(i)*c(i-1)/d1(i-1);
            b1(i)=b(i)-a(i)*b1(i-1)/d1(i-1);
        else
            a(i)=1/(x(i)-x(i-1));
            d(i)=-(1/(x(i+1)-x(i))+1/(x(i)-x(i-1)));
            c(i)=1/(x(i+1)-x(i));
            b(i)=lap_phi(i)*cell_vol(i);
            %
            d1(i)=d(i)-a(i)*c(i-1)/d1(i-1);
            b1(i)=b(i)-a(i)*b1(i-1)/d1(i-1);
            A(i,i)=d(i);
        end
    end
    % Algoritmo Thomas - Substituição Ascendente %
    for i=cell_num:-1:1
        if i==cell_num
            phi_num(i)=b1(i)/d1(i);
        else
            phi_num(i)=(b1(i)-c(i)*phi_num(i+1))/d1(i);
        end
    end
    TempoInverterMatriz=cputime-TempoInverterMatriz;
    
    
    % Cálculo do Laplaciano Numerico %
    for i=1:cell_num
        if i==1
            if dirichlet
                flux_phi_west=(phi_num(i)-phi_West)/x(i);
                flux_phi_east=(phi_num(i+1)-phi_num(i))/(x(i+1)-x(i));
            else
                flux_phi_west=flux_phi_West;
                flux_phi_east=(phi_num(i+1)-phi_num(i))/(x(i+1)-x(i));
            end
        elseif i==cell_num
            flux_phi_west=(phi_num(i)-phi_num(i-1))/(x(i)-x(i-1));
            flux_phi_east=(phi_East-phi_num(i-1))/(L-x(i));
        else 
            flux_phi_west=(phi_num(i)-phi_num(i-1))/(x(i)-x(i-1));
            flux_phi_east=(phi_num(i+1)-phi_num(i))/(x(i+1)-x(i));
        end
        lap_phi_num(i)=(flux_phi_east-flux_phi_west)/cell_vol(i);
    end
end
% Fim da Função %
end
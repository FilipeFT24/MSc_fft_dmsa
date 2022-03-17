function [phi_num,lap_phi_num]=Scheme6thOrder(explicito,dirichlet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   16 de Junho de 2016                                             %
%                                   28 de Junho de 2016                                             %
%                                                                                                   %
% Resolução do Problema de 6ª Ordem                                                                 %
%                                                                                                   %
% Condições de Fronteira de Dirichlet nas duas fronteiras;                                          %
% Condições de Fronteira de Neumann na fronteira esquerda e de Dirichlet na fronteira esquerda;     %
%                                                                                                   %
% Cálculo Explicito é a verificação dos fluxos, reconstroi-se o polinomio nas faces com base em     %
% diferenças finitas de 6ª ordem e depois a partir do método de volume finito calcula-se os fluxos  %
% nas faces da célula.                                                                              %
%                                                                                                   %
% Cálculo Implicito é a determinação da propriedade no centroide da célula com base nas             %
% reconstruçoes realizadas para as faces, sendo utilizado o método de volume finito.                %
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
%% Definição dos Fluxos nas Faces %%
%
% Cálculo das Distancias entre a face e o centroide das celulas utilizadas no stencil %
%
for i=1:face_num
    if i==1
        d1=0;
        d2=x(i)-faces(i);
        d3=x(i+1)-faces(i);
        d4=x(i+2)-faces(i);
        d5=x(i+3)-faces(i);
        d6=x(i+4)-faces(i);
    elseif i==2        
        d1=x(i)-faces(i);
        d2=x(i+1)-faces(i);
        d3=x(i+2)-faces(i);
        d4=x(i+3)-faces(i);
        d5=x(i-1)-faces(i);
        d6=faces(1)-faces(i);
    elseif i==3
        d1=x(i)-faces(i);
        d2=x(i+1)-faces(i);
        d3=x(i+2)-faces(i);
        d4=x(i-1)-faces(i);
        d5=x(i-2)-faces(i);
        d6=faces(1)-faces(i);       
    elseif i==face_num-2
        d1=x(i)-faces(i);
        d2=x(i+1)-faces(i);
        d3=faces(face_num)-faces(i);
        d4=x(i-1)-faces(i);
        d5=x(i-2)-faces(i);
        d6=x(i-3)-faces(i);
    elseif i==face_num-1
        d1=x(i)-faces(i);
        d2=faces(face_num)-faces(i);
        d3=x(i-1)-faces(i);
        d4=x(i-2)-faces(i);
        d5=x(i-3)-faces(i);
        d6=x(i-4)-faces(i);
    elseif i==face_num
        d1=0;
        d2=x(i-1)-faces(i);
        d3=x(i-2)-faces(i);
        d4=x(i-3)-faces(i);
        d5=x(i-4)-faces(i);
        d6=x(i-5)-faces(i);
    else
        d1=x(i)-faces(i);
        d2=x(i+1)-faces(i);
        d3=x(i+2)-faces(i);
        d4=x(i-1)-faces(i);
        d5=x(i-2)-faces(i);
        d6=x(i-3)-faces(i);
    end
    %
    % Inversão da Matriz Para se obter os Coeficientes da expressão de diferenças finitas %
    %
    if i==2 && ~ dirichlet || i==3 && ~ dirichlet
        C=[ 1        1        1        1        1        0      ;
            d1       d2       d3       d4       d5       1      ;
            d1^2/2   d2^2/2   d3^2/2   d4^2/2   d5^2/2   d6^1   ;
            d1^3/6   d2^3/6   d3^3/6   d4^3/6   d5^3/6   d6^2/2 ;
            d1^4/24  d2^4/24  d3^4/24  d4^4/24  d5^4/24  d6^3/6 ;
            d1^5/120 d2^5/120 d3^5/120 d4^5/120 d5^5/120 d6^4/24];        
    else
        C=[ 1    1    1    1    1    1  ;
            d1   d2   d3   d4   d5   d6  ;
            d1^2 d2^2 d3^2 d4^2 d5^2 d6^2;
            d1^3 d2^3 d3^3 d4^3 d5^3 d6^3;
            d1^4 d2^4 d3^4 d4^4 d5^4 d6^4;
            d1^5 d2^5 d3^5 d4^5 d5^5 d6^5];
    end
    e=[0;1;0;0;0;0];
    %
    const=inv(C)*e;
    %
    % Valores dos Coeficientes do FDM %
    %
    coefs(i,1)=const(1);
    coefs(i,2)=const(2);
    coefs(i,3)=const(3);
    coefs(i,4)=const(4);
    coefs(i,5)=const(5);
    coefs(i,6)=const(6);
end
%
%% Metodo Explicito %%
%
if explicito
    TempoInverterMatriz=cputime;
    for i=1:cell_num
        if i==1
            if dirichlet
                flux_phi_east=coefs(i+1,1)*phi(i+1)+coefs(i+1,2)*phi(i+2)+coefs(i+1,3)*phi(i+3)+coefs(i+1,4)*phi(i+4)+coefs(i+1,5)*phi(i+0)+coefs(i+1,6)*phi_West;
                flux_phi_west=coefs(i+0,1)*phi_West+coefs(i+0,2)*phi(i+0)+coefs(i+0,3)*phi(i+1)+coefs(i+0,4)*phi(i+2)+coefs(i+0,5)*phi(i+3)+coefs(i+0,6)*phi(i+4);
            else
                flux_phi_east=coefs(i+1,1)*phi(i+1)+coefs(i+1,2)*phi(i+2)+coefs(i+1,3)*phi(i+3)+coefs(i+1,4)*phi(i+4)+coefs(i+1,5)*phi(i+0)+coefs(i+1,6)*flux_phi_West;
                flux_phi_west=flux_phi_West;
            end
        elseif i==2
            if dirichlet
                flux_phi_east=coefs(i+1,1)*phi(i+1)+coefs(i+1,2)*phi(i+2)+coefs(i+1,3)*phi(i+3)+coefs(i+1,4)*phi(i+0)+coefs(i+1,5)*phi(i-1)+coefs(i+1,6)*phi_West;
                flux_phi_west=coefs(i+0,1)*phi(i+0)+coefs(i+0,2)*phi(i+1)+coefs(i+0,3)*phi(i+2)+coefs(i+0,4)*phi(i+3)+coefs(i+0,5)*phi(i-1)+coefs(i+0,6)*phi_West;
            else
                flux_phi_east=coefs(i+1,1)*phi(i+1)+coefs(i+1,2)*phi(i+2)+coefs(i+1,3)*phi(i+3)+coefs(i+1,4)*phi(i+0)+coefs(i+1,5)*phi(i-1)+coefs(i+1,6)*flux_phi_West;
                flux_phi_west=coefs(i+0,1)*phi(i+0)+coefs(i+0,2)*phi(i+1)+coefs(i+0,3)*phi(i+2)+coefs(i+0,4)*phi(i+3)+coefs(i+0,5)*phi(i-1)+coefs(i+0,6)*flux_phi_West;
            end
        elseif i==3
            if dirichlet
                flux_phi_east=coefs(i+1,1)*phi(i+1)+coefs(i+1,2)*phi(i+2)+coefs(i+1,3)*phi(i+3)+coefs(i+1,4)*phi(i+0)+coefs(i+1,5)*phi(i-1)+coefs(i+1,6)*phi(i-2);
                flux_phi_west=coefs(i+0,1)*phi(i+0)+coefs(i+0,2)*phi(i+1)+coefs(i+0,3)*phi(i+2)+coefs(i+0,4)*phi(i-1)+coefs(i+0,5)*phi(i-2)+coefs(i+0,6)*phi_West;
            else
                flux_phi_east=coefs(i+1,1)*phi(i+1)+coefs(i+1,2)*phi(i+2)+coefs(i+1,3)*phi(i+3)+coefs(i+1,4)*phi(i+0)+coefs(i+1,5)*phi(i-1)+coefs(i+1,6)*phi(i-2);
                flux_phi_west=coefs(i+0,1)*phi(i+0)+coefs(i+0,2)*phi(i+1)+coefs(i+0,3)*phi(i+2)+coefs(i+0,4)*phi(i-1)+coefs(i+0,5)*phi(i-2)+coefs(i+0,6)*flux_phi_West;
            end
        elseif i==cell_num-2
            flux_phi_east=coefs(i+1,1)*phi(i+1)+coefs(i+1,2)*phi(i+2)+coefs(i+1,3)*phi_East+coefs(i+1,4)*phi(i+0)+coefs(i+1,5)*phi(i-1)+coefs(i+1,6)*phi(i-2);
            flux_phi_west=coefs(i+0,1)*phi(i+0)+coefs(i+0,2)*phi(i+1)+coefs(i+0,3)*phi(i+2)+coefs(i+0,4)*phi(i-1)+coefs(i+0,5)*phi(i-2)+coefs(i+0,6)*phi(i-3);
        elseif i==cell_num-1
            flux_phi_east=coefs(i+1,1)*phi(i+1)+coefs(i+1,2)*phi_East+coefs(i+1,3)*phi(i+0)+coefs(i+1,4)*phi(i-1)+coefs(i+1,5)*phi(i-2)+coefs(i+1,6)*phi(i-3);
            flux_phi_west=coefs(i+0,1)*phi(i+0)+coefs(i+0,2)*phi(i+1)+coefs(i+0,3)*phi_East+coefs(i+0,4)*phi(i-1)+coefs(i+0,5)*phi(i-2)+coefs(i+0,6)*phi(i-3);
        elseif i==cell_num
            flux_phi_east=coefs(i+1,1)*phi_East+coefs(i+1,2)*phi(i+0)+coefs(i+1,3)*phi(i-1)+coefs(i+1,4)*phi(i-2)+coefs(i+1,5)*phi(i-3)+coefs(i+1,6)*phi(i-4);
            flux_phi_west=coefs(i+0,1)*phi(i+0)+coefs(i+0,2)*phi_East+coefs(i+0,3)*phi(i-1)+coefs(i+0,4)*phi(i-2)+coefs(i+0,5)*phi(i-3)+coefs(i+0,6)*phi(i-4);
        else
            flux_phi_east=coefs(i+1,1)*phi(i+1)+coefs(i+1,2)*phi(i+2)+coefs(i+1,3)*phi(i+3)+coefs(i+1,4)*phi(i+0)+coefs(i+1,5)*phi(i-1)+coefs(i+1,6)*phi(i-2);
            flux_phi_west=coefs(i+0,1)*phi(i+0)+coefs(i+0,2)*phi(i+1)+coefs(i+0,3)*phi(i+2)+coefs(i+0,4)*phi(i-1)+coefs(i+0,5)*phi(i-2)+coefs(i+0,6)*phi(i-3);
        end
        lap_phi_num(i)=(flux_phi_east-flux_phi_west)/cell_vol(i);
        phi_num(i)=phi(i);
    end
    TempoInverterMatriz=cputime-TempoInverterMatriz;
    %% Mátodo Implicito %%
else
    TempoInverterMatriz=cputime;
    A=zeros(cell_num);
    for i=1:cell_num
        if i==1
            if dirichlet
                A(i,i)=coefs(i+1,5)-coefs(i+0,2);
                A(i,i+1)=coefs(i+1,1)-coefs(i+0,3);
                A(i,i+2)=coefs(i+1,2)-coefs(i+0,4);
                A(i,i+3)=coefs(i+1,3)-coefs(i+0,5);
                A(i,i+4)=coefs(i+1,4)-coefs(i+0,6);
                b(i)=lap_phi(i)*cell_vol(i)-(coefs(i+1,6)-coefs(i+0,1))*phi_West;
            else
                A(i,i)=coefs(i+1,5);
                A(i,i+1)=coefs(i+1,1);
                A(i,i+2)=coefs(i+1,2);
                A(i,i+3)=coefs(i+1,3);
                A(i,i+4)=coefs(i+1,4);
                b(i)=lap_phi(i)*cell_vol(i)-(coefs(i+1,6)-1)*flux_phi_West;
            end
        elseif i==2
            if dirichlet
                A(i,i-1)=coefs(i+1,5)-coefs(i+0,5);
                A(i,i)=coefs(i+1,4)-coefs(i+0,1);
                A(i,i+1)=coefs(i+1,1)-coefs(i+0,2);
                A(i,i+2)=coefs(i+1,2)-coefs(i+0,3);
                A(i,i+3)=coefs(i+1,3)-coefs(i+0,4);
                b(i)=lap_phi(i)*cell_vol(i)-(coefs(i+1,6)-coefs(i+0,6))*phi_West;
            else
                A(i,i-1)=coefs(i+1,5)-coefs(i+0,5);
                A(i,i)=coefs(i+1,4)-coefs(i+0,1);
                A(i,i+1)=coefs(i+1,1)-coefs(i+0,2);
                A(i,i+2)=coefs(i+1,2)-coefs(i+0,3);
                A(i,i+3)=coefs(i+1,3)-coefs(i+0,4);
                b(i)=lap_phi(i)*cell_vol(i)-(coefs(i+1,6)-coefs(i+0,6))*flux_phi_West;
            end
        elseif i==3
            if dirichlet
                A(i,i-2)=coefs(i+1,6)-coefs(i+0,5);
                A(i,i-1)=coefs(i+1,5)-coefs(i+0,4);
                A(i,i)=coefs(i+1,4)-coefs(i+0,1);
                A(i,i+1)=coefs(i+1,1)-coefs(i+0,2);
                A(i,i+2)=coefs(i+1,2)-coefs(i+0,3);
                A(i,i+3)=coefs(i+1,3);
                b(i)=lap_phi(i)*cell_vol(i)+coefs(i+0,6)*phi_West;
            else
                A(i,i-2)=coefs(i+1,6)-coefs(i+0,5);
                A(i,i-1)=coefs(i+1,5)-coefs(i+0,4);
                A(i,i)=coefs(i+1,4)-coefs(i+0,1);
                A(i,i+1)=coefs(i+1,1)-coefs(i+0,2);
                A(i,i+2)=coefs(i+1,2)-coefs(i+0,3);
                A(i,i+3)=coefs(i+1,3);
                b(i)=lap_phi(i)*cell_vol(i)+coefs(i+0,6)*flux_phi_West;
            end
        elseif i==cell_num-2
            A(i,i-3)=-coefs(i+0,6);
            A(i,i-2)=coefs(i+1,6)-coefs(i+0,5);
            A(i,i-1)=coefs(i+1,5)-coefs(i+0,4);
            A(i,i)=coefs(i+1,4)-coefs(i+0,1);
            A(i,i+1)=coefs(i+1,1)-coefs(i+0,2);
            A(i,i+2)=coefs(i+1,2)-coefs(i+0,3);
            b(i)=lap_phi(i)*cell_vol(i)-coefs(i+1,3)*phi_East;
        elseif i==cell_num-1
            A(i,i-3)=coefs(i+1,6)-coefs(i+0,6);
            A(i,i-2)=coefs(i+1,5)-coefs(i+0,5);
            A(i,i-1)=coefs(i+1,4)-coefs(i+0,4);
            A(i,i)=coefs(i+1,3)-coefs(i+0,1);
            A(i,i+1)=coefs(i+1,1)-coefs(i+0,2);
            b(i)=lap_phi(i)*cell_vol(i)-(coefs(i+1,2)-coefs(i+0,3))*phi_East;
        elseif i==cell_num
            A(i,i-4)=coefs(i+1,6)-coefs(i+0,6);
            A(i,i-3)=coefs(i+1,5)-coefs(i+0,5);
            A(i,i-2)=coefs(i+1,4)-coefs(i+0,4);
            A(i,i-1)=coefs(i+1,3)-coefs(i+0,3);
            A(i,i)=coefs(i+1,2)-coefs(i+0,1);
            b(i)=lap_phi(i)*cell_vol(i)-(coefs(i+1,1)-coefs(i+0,2))*phi_East;
        else
            A(i,i-3)=-coefs(i+0,6);
            A(i,i-2)=coefs(i+1,6)-coefs(i+0,5);
            A(i,i-1)=coefs(i+1,5)-coefs(i+0,4);
            A(i,i)=coefs(i+1,4)-coefs(i+0,1);
            A(i,i+1)=coefs(i+1,1)-coefs(i+0,2);
            A(i,i+2)=coefs(i+1,2)-coefs(i+0,3);
            A(i,i+3)=coefs(i+1,3);
            b(i)=lap_phi(i)*cell_vol(i);
        end
    end
    phi_num=inv(A)*b';
    TempoInverterMatriz=cputime-TempoInverterMatriz;
    for i=1:cell_num
        if i==1
            if dirichlet
                flux_phi_east=coefs(i+1,1)*phi_num(i+1)+coefs(i+1,2)*phi_num(i+2)+coefs(i+1,3)*phi_num(i+3)+coefs(i+1,4)*phi_num(i+4)+coefs(i+1,5)*phi_num(i+0)+coefs(i+1,6)*phi_West;
                flux_phi_west=coefs(i+0,1)*phi_West+coefs(i+0,2)*phi_num(i+0)+coefs(i+0,3)*phi_num(i+1)+coefs(i+0,4)*phi_num(i+2)+coefs(i+0,5)*phi_num(i+3)+coefs(i+0,6)*phi_num(i+4);
            else
                flux_phi_east=coefs(i+1,1)*phi_num(i+1)+coefs(i+1,2)*phi_num(i+2)+coefs(i+1,3)*phi_num(i+3)+coefs(i+1,4)*phi_num(i+4)+coefs(i+1,5)*phi_num(i+0)+coefs(i+1,6)*flux_phi_West;
                flux_phi_west=flux_phi_West;
            end
        elseif i==2
            if dirichlet
                flux_phi_east=coefs(i+1,1)*phi_num(i+1)+coefs(i+1,2)*phi_num(i+2)+coefs(i+1,3)*phi_num(i+3)+coefs(i+1,4)*phi_num(i+0)+coefs(i+1,5)*phi_num(i-1)+coefs(i+1,6)*phi_West;
                flux_phi_west=coefs(i+0,1)*phi_num(i+0)+coefs(i+0,2)*phi_num(i+1)+coefs(i+0,3)*phi_num(i+2)+coefs(i+0,4)*phi_num(i+3)+coefs(i+0,5)*phi_num(i-1)+coefs(i+0,6)*phi_West;
            else
                flux_phi_east=coefs(i+1,1)*phi_num(i+1)+coefs(i+1,2)*phi_num(i+2)+coefs(i+1,3)*phi_num(i+3)+coefs(i+1,4)*phi_num(i+0)+coefs(i+1,5)*phi_num(i-1)+coefs(i+1,6)*flux_phi_West;
                flux_phi_west=coefs(i+0,1)*phi_num(i+0)+coefs(i+0,2)*phi_num(i+1)+coefs(i+0,3)*phi_num(i+2)+coefs(i+0,4)*phi_num(i+3)+coefs(i+0,5)*phi_num(i-1)+coefs(i+0,6)*flux_phi_West;
            end
        elseif i==3
            if dirichlet
                flux_phi_east=coefs(i+1,1)*phi_num(i+1)+coefs(i+1,2)*phi_num(i+2)+coefs(i+1,3)*phi_num(i+3)+coefs(i+1,4)*phi_num(i+0)+coefs(i+1,5)*phi_num(i-1)+coefs(i+1,6)*phi_num(i-2);
                flux_phi_west=coefs(i+0,1)*phi_num(i+0)+coefs(i+0,2)*phi_num(i+1)+coefs(i+0,3)*phi_num(i+2)+coefs(i+0,4)*phi_num(i-1)+coefs(i+0,5)*phi_num(i-2)+coefs(i+0,6)*phi_West;
            else
                flux_phi_east=coefs(i+1,1)*phi_num(i+1)+coefs(i+1,2)*phi_num(i+2)+coefs(i+1,3)*phi_num(i+3)+coefs(i+1,4)*phi_num(i+0)+coefs(i+1,5)*phi_num(i-1)+coefs(i+1,6)*phi_num(i-2);
                flux_phi_west=coefs(i+0,1)*phi_num(i+0)+coefs(i+0,2)*phi_num(i+1)+coefs(i+0,3)*phi_num(i+2)+coefs(i+0,4)*phi_num(i-1)+coefs(i+0,5)*phi_num(i-2)+coefs(i+0,6)*flux_phi_West;
            end
        elseif i==cell_num-2
            flux_phi_east=coefs(i+1,1)*phi_num(i+1)+coefs(i+1,2)*phi_num(i+2)+coefs(i+1,3)*phi_East+coefs(i+1,4)*phi_num(i+0)+coefs(i+1,5)*phi_num(i-1)+coefs(i+1,6)*phi_num(i-2);
            flux_phi_west=coefs(i+0,1)*phi_num(i+0)+coefs(i+0,2)*phi_num(i+1)+coefs(i+0,3)*phi_num(i+2)+coefs(i+0,4)*phi_num(i-1)+coefs(i+0,5)*phi_num(i-2)+coefs(i+0,6)*phi_num(i-3);
        elseif i==cell_num-1
            flux_phi_east=coefs(i+1,1)*phi_num(i+1)+coefs(i+1,2)*phi_East+coefs(i+1,3)*phi_num(i+0)+coefs(i+1,4)*phi_num(i-1)+coefs(i+1,5)*phi_num(i-2)+coefs(i+1,6)*phi_num(i-3);
            flux_phi_west=coefs(i+0,1)*phi_num(i+0)+coefs(i+0,2)*phi_num(i+1)+coefs(i+0,3)*phi_East+coefs(i+0,4)*phi_num(i-1)+coefs(i+0,5)*phi_num(i-2)+coefs(i+0,6)*phi_num(i-3);
        elseif i==cell_num
            flux_phi_east=coefs(i+1,1)*phi_East+coefs(i+1,2)*phi_num(i+0)+coefs(i+1,3)*phi_num(i-1)+coefs(i+1,4)*phi_num(i-2)+coefs(i+1,5)*phi_num(i-3)+coefs(i+1,6)*phi_num(i-4);
            flux_phi_west=coefs(i+0,1)*phi_num(i+0)+coefs(i+0,2)*phi_East+coefs(i+0,3)*phi_num(i-1)+coefs(i+0,4)*phi_num(i-2)+coefs(i+0,5)*phi_num(i-3)+coefs(i+0,6)*phi_num(i-4);
        else
            flux_phi_east=coefs(i+1,1)*phi_num(i+1)+coefs(i+1,2)*phi_num(i+2)+coefs(i+1,3)*phi_num(i+3)+coefs(i+1,4)*phi_num(i+0)+coefs(i+1,5)*phi_num(i-1)+coefs(i+1,6)*phi_num(i-2);
            flux_phi_west=coefs(i+0,1)*phi_num(i+0)+coefs(i+0,2)*phi_num(i+1)+coefs(i+0,3)*phi_num(i+2)+coefs(i+0,4)*phi_num(i-1)+coefs(i+0,5)*phi_num(i-2)+coefs(i+0,6)*phi_num(i-3);
        end
        lap_phi_num(i)=(flux_phi_east-flux_phi_west)/cell_vol(i);
    end
end

end

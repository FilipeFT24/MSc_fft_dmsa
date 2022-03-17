function [phi_num,lap_phi_num]=Scheme4thOrder(explicito,dirichlet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   08 de Junho de 2016                                             %
%                                   28 de Junho de 2016                                             %
%                                                                                                   %
% Resolu��o do Problema de 4� Ordem                                                                 %
%                                                                                                   %
% Condi��es de Fronteira de Dirichlet nas duas fronteiras;                                          %
% Condi��es de Fronteira de Neumann na fronteira esquerda e de Dirichlet na fronteira esquerda;     %
%                                                                                                   %
% C�lculo Explicito � a verifica��o dos fluxos, reconstroi-se o polinomio nas faces com base em     %
% diferen�as finitas de 4� ordem e depois a partir do m�todo de volume finito calcula-se os fluxos  %
% nas faces da c�lula.                                                                              %
%                                                                                                   %
% C�lculo Implicito � a determina��o da propriedade no centroide da c�lula com base nas             %
% reconstru�oes realizadas para as faces, sendo utilizado o m�todo de volume finito.                %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Declara��o das Variaveis Globais %%
%
global TempoInverterMatriz;
global L Lref cell_num face_num;
global x faces cell_vol;
global phi lap_phi phi_West phi_East flux_phi_West flux_phi_East;
%
%% Defini��o dos Fluxos nas Faces %%
%
for i=1:face_num
    if i==1
        d1=0;
        d2=x(i)-faces(i);
        d3=x(i+1)-faces(i);
        d4=x(i+2)-faces(i);
    elseif i==2
        d1=x(i)-faces(i);
        d2=x(i+1)-faces(i);
        d3=x(i-1)-faces(i);
        d4=faces(1)-faces(i);
    elseif i==face_num-1
        d1=x(i)-faces(i);
        d2=faces(face_num)-faces(i);
        d3=x(i-1)-faces(i);
        d4=x(i-2)-faces(i);
    elseif i==face_num
        d1=0;
        d2=x(i-1)-faces(i);
        d3=x(i-2)-faces(i);
        d4=x(i-3)-faces(i);
    else
        d1=x(i)-faces(i);
        d2=x(i+1)-faces(i);
        d3=x(i-1)-faces(i);
        d4=x(i-2)-faces(i);
    end
    %
    % Valores dos Coeficientes do FDM para as Condi��es de Dirichlet %
    %
    if i==1 && ~ dirichlet
        coefs(i,1)=0;
        coefs(i,2)=0;
        coefs(i,3)=0;
        coefs(i,4)=0;
    elseif i==2 && ~ dirichlet
        coefs(i,1)=(2*d2^2*d4 + 2*d2*d3*d4 - 3*d2*d4^2 + 2*d3^2*d4 - 3*d3*d4^2)/(d1^3*d2 + d1^3*d3 - 2*d1^3*d4 - d1^2*d2^2 - d1^2*d2*d3 - d1^2*d3^2 + 3*d1^2*d4^2 + 2*d1*d2^2*d4 + 2*d1*d2*d3*d4 - 3*d1*d2*d4^2 + 2*d1*d3^2*d4 - 3*d1*d3*d4^2 + d2^2*d3^2 - 2*d2^2*d3*d4 - 2*d2*d3^2*d4 + 3*d2*d3*d4^2);
        coefs(i,2)=(2*d1^2*d4 + 2*d1*d3*d4 - 3*d1*d4^2 + 2*d3^2*d4 - 3*d3*d4^2)/(- d1^2*d2^2 + 2*d1^2*d2*d4 + d1^2*d3^2 - 2*d1^2*d3*d4 + d1*d2^3 - d1*d2^2*d3 + 2*d1*d2*d3*d4 - 3*d1*d2*d4^2 - 2*d1*d3^2*d4 + 3*d1*d3*d4^2 + d2^3*d3 - 2*d2^3*d4 - d2^2*d3^2 + 3*d2^2*d4^2 + 2*d2*d3^2*d4 - 3*d2*d3*d4^2);
        coefs(i,3)=(2*d1^2*d4 + 2*d1*d2*d4 - 3*d1*d4^2 + 2*d2^2*d4 - 3*d2*d4^2)/(d1^2*d2^2 - 2*d1^2*d2*d4 - d1^2*d3^2 + 2*d1^2*d3*d4 - 2*d1*d2^2*d4 - d1*d2*d3^2 + 2*d1*d2*d3*d4 + 3*d1*d2*d4^2 + d1*d3^3 - 3*d1*d3*d4^2 - d2^2*d3^2 + 2*d2^2*d3*d4 + d2*d3^3 - 3*d2*d3*d4^2 - 2*d3^3*d4 + 3*d3^2*d4^2);
        coefs(i,4)=(d1*d2 + d1*d3 + d2*d3)/(d1*d2 + d1*d3 - 2*d1*d4 + d2*d3 - 2*d2*d4 - 2*d3*d4 + 3*d4^2);
    else
        coefs(i,1)=-(d2*d3 + d2*d4 + d3*d4)/(d1^2*d2 + d1^2*d3 + d1^2*d4 - d1^3 - d1*d2*d3 - d1*d2*d4 - d1*d3*d4 + d2*d3*d4);
        coefs(i,2)=-(d1*d3 + d1*d4 + d3*d4)/(d1*d2^2 + d2^2*d3 + d2^2*d4 - d2^3 - d1*d2*d3 - d1*d2*d4 + d1*d3*d4 - d2*d3*d4);
        coefs(i,3)=-(d1*d2 + d1*d4 + d2*d4)/(d1*d3^2 + d2*d3^2 + d3^2*d4 - d3^3 - d1*d2*d3 + d1*d2*d4 - d1*d3*d4 - d2*d3*d4);
        coefs(i,4)=-(d1*d2 + d1*d3 + d2*d3)/(d1*d4^2 + d2*d4^2 + d3*d4^2 - d4^3 + d1*d2*d3 - d1*d2*d4 - d1*d3*d4 - d2*d3*d4);
    end
end
%
%% Metodo %%
%
if explicito
    TempoInverterMatriz=cputime;
    for i=1:cell_num
        if i==1
            if dirichlet
                flux_phi_east=coefs(i+1,1)*phi(i+1)+coefs(i+1,2)*phi(i+2)+coefs(i+1,3)*phi(i)+coefs(i+1,4)*phi_West;
                flux_phi_west=coefs(i,1)*phi_West+coefs(i,2)*phi(i)+coefs(i,3)*phi(i+1)+coefs(i,4)*phi(i+2);
            else
                flux_phi_east=coefs(i+1,1)*phi(i+1)+coefs(i+1,2)*phi(i+2)+coefs(i+1,3)*phi(i)+coefs(i+1,4)*flux_phi_West;
                flux_phi_west=flux_phi_West;
            end
        elseif i==2
            if dirichlet
                flux_phi_east=coefs(i+1,1)*phi(i+1)+coefs(i+1,2)*phi(i+2)+coefs(i+1,3)*phi(i)+coefs(i+1,4)*phi(i-1);
                flux_phi_west=coefs(i,1)*phi(i)+coefs(i,2)*phi(i+1)+coefs(i,3)*phi(i-1)+coefs(i,4)*phi_West;
            else
                flux_phi_east=coefs(i+1,1)*phi(i+1)+coefs(i+1,2)*phi(i+2)+coefs(i+1,3)*phi(i)+coefs(i+1,4)*phi(i-1);
                flux_phi_west=coefs(i,1)*phi(i)+coefs(i,2)*phi(i+1)+coefs(i,3)*phi(i-1)+coefs(i,4)*flux_phi_West;
            end
        elseif i==cell_num-1
            flux_phi_east=coefs(i+1,1)*phi(i+1)+coefs(i+1,2)*phi_East+coefs(i+1,3)*phi(i)+coefs(i+1,4)*phi(i-1);
            flux_phi_west=coefs(i,1)*phi(i)+coefs(i,2)*phi(i+1)+coefs(i,3)*phi(i-1)+coefs(i,4)*phi(i-2);
        elseif i==cell_num
            flux_phi_east=coefs(i+1,1)*phi_East+coefs(i+1,2)*phi(i)+coefs(i+1,3)*phi(i-1)+coefs(i+1,4)*phi(i-2);
            flux_phi_west=coefs(i,1)*phi(i)+coefs(i,2)*phi_East+coefs(i,3)*phi(i-1)+coefs(i,4)*phi(i-2);
        else
            flux_phi_east=coefs(i+1,1)*phi(i+1)+coefs(i+1,2)*phi(i+2)+coefs(i+1,3)*phi(i)+coefs(i+1,4)*phi(i-1);
            flux_phi_west=coefs(i,1)*phi(i)+coefs(i,2)*phi(i+1)+coefs(i,3)*phi(i-1)+coefs(i,4)*phi(i-2);
        end
        lap_phi_num(i)=(flux_phi_east-flux_phi_west)/cell_vol(i);
        phi_num(i)=phi(i);
    end
    TempoInverterMatriz=cputime-TempoInverterMatriz;
else
    TempoInverterMatriz=cputime;
    A=zeros(cell_num);
    for i=1:cell_num
        if i==1
            if dirichlet
                A(i,i)=coefs(i+1,3)-coefs(i,2);
                A(i,i+1)=coefs(i+1,1)-coefs(i,3);
                A(i,i+2)=coefs(i+1,2)-coefs(i,4);
                b(i)=lap_phi(i)*cell_vol(i)+(coefs(i,1)-coefs(i+1,4))*phi_West;
            else
                A(i,i)=coefs(i+1,3);
                A(i,i+1)=coefs(i+1,1);
                A(i,i+2)=coefs(i+1,2);
                b(i)=lap_phi(i)*cell_vol(i)+(1-coefs(i+1,4))*flux_phi_West;
            end
        elseif i==2
            if dirichlet
                A(i,i-1)=coefs(i+1,4)-coefs(i,3);
                A(i,i)=coefs(i+1,3)-coefs(i,1);
                A(i,i+1)=coefs(i+1,1)-coefs(i,2);
                A(i,i+2)=coefs(i+1,2);
                b(i)=lap_phi(i)*cell_vol(i)+coefs(i,4)*phi_West;
            else
                A(i,i-1)=coefs(i+1,4)-coefs(i,3);
                A(i,i)=coefs(i+1,3)-coefs(i,1);
                A(i,i+1)=coefs(i+1,1)-coefs(i,2);
                A(i,i+2)=coefs(i+1,2);
                b(i)=lap_phi(i)*cell_vol(i)+coefs(i,4)*flux_phi_West;
            end
        elseif i==cell_num-1
            A(i,i-2)=-coefs(i,4);
            A(i,i-1)=coefs(i+1,4)-coefs(i,3);
            A(i,i)=coefs(i+1,3)-coefs(i,1);
            A(i,i+1)=coefs(i+1,1)-coefs(i,2);
            b(i)=lap_phi(i)*cell_vol(i)-coefs(i+1,2)*phi_East;
        elseif i==cell_num
            A(i,i-2)=coefs(i+1,4)-coefs(i,4);
            A(i,i-1)=coefs(i+1,3)-coefs(i,3);
            A(i,i)=coefs(i+1,2)-coefs(i,1);
            b(i)=lap_phi(i)*cell_vol(i)+(coefs(i,2)-coefs(i+1,1))*phi_East;
        else
            A(i,i-2)=-coefs(i,4);
            A(i,i-1)=coefs(i+1,4)-coefs(i,3);
            A(i,i)=coefs(i+1,3)-coefs(i,1);
            A(i,i+1)=coefs(i+1,1)-coefs(i,2);
            A(i,i+2)=coefs(i+1,2);
            b(i)=lap_phi(i)*cell_vol(i);
        end
    end
    phi_num=inv(A)*b';
    TempoInverterMatriz=cputime-TempoInverterMatriz;
    for i=1:cell_num
        if i==1
            if dirichlet
                flux_phi_east=coefs(i+1,1)*phi_num(i+1)+coefs(i+1,2)*phi_num(i+2)+coefs(i+1,3)*phi_num(i)+coefs(i+1,4)*phi_West;
                flux_phi_west=coefs(i,1)*phi_West+coefs(i,2)*phi_num(i)+coefs(i,3)*phi_num(i+1)+coefs(i,4)*phi_num(i+2);
            else
                flux_phi_east=coefs(i+1,1)*phi_num(i+1)+coefs(i+1,2)*phi_num(i+2)+coefs(i+1,3)*phi_num(i)+coefs(i+1,4)*flux_phi_West;
                flux_phi_west=flux_phi_West;
            end
        elseif i==2
            if dirichlet
                flux_phi_east=coefs(i+1,1)*phi_num(i+1)+coefs(i+1,2)*phi_num(i+2)+coefs(i+1,3)*phi_num(i)+coefs(i+1,4)*phi_num(i-1);
                flux_phi_west=coefs(i,1)*phi_num(i)+coefs(i,2)*phi_num(i+1)+coefs(i,3)*phi_num(i-1)+coefs(i,4)*phi_West;
            else
                flux_phi_east=coefs(i+1,1)*phi_num(i+1)+coefs(i+1,2)*phi_num(i+2)+coefs(i+1,3)*phi_num(i)+coefs(i+1,4)*phi_num(i-1);
                flux_phi_west=coefs(i,1)*phi_num(i)+coefs(i,2)*phi_num(i+1)+coefs(i,3)*phi_num(i-1)+coefs(i,4)*flux_phi_West;
            end
        elseif i==cell_num-1
            flux_phi_east=coefs(i+1,1)*phi_num(i+1)+coefs(i+1,2)*phi_East+coefs(i+1,3)*phi_num(i)+coefs(i+1,4)*phi_num(i-1);
            flux_phi_west=coefs(i,1)*phi_num(i)+coefs(i,2)*phi_num(i+1)+coefs(i,3)*phi_num(i-1)+coefs(i,4)*phi_num(i-2);
        elseif i==cell_num
            flux_phi_east=coefs(i+1,1)*phi_East+coefs(i+1,2)*phi_num(i)+coefs(i+1,3)*phi_num(i-1)+coefs(i+1,4)*phi_num(i-2);
            flux_phi_west=coefs(i,1)*phi_num(i)+coefs(i,2)*phi_East+coefs(i,3)*phi_num(i-1)+coefs(i,4)*phi_num(i-2);
        else
            flux_phi_east=coefs(i+1,1)*phi_num(i+1)+coefs(i+1,2)*phi_num(i+2)+coefs(i+1,3)*phi_num(i)+coefs(i+1,4)*phi_num(i-1);
            flux_phi_west=coefs(i,1)*phi_num(i)+coefs(i,2)*phi_num(i+1)+coefs(i,3)*phi_num(i-1)+coefs(i,4)*phi_num(i-2);
        end
        lap_phi_num(i)=(flux_phi_east-flux_phi_west)/cell_vol(i);
    end
end

end
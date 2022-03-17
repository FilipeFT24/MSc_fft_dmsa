function [A,b]=matrix2ndOrder_convection
%
%% Declaração das Variaveis Globais %%
%
global L Lref cell_num face_num;
global x faces cell_vol;
global phi lap_phi phi_West phi_East;
%
%% Construção da Matriz %%
%
% Cálculo das Distancias entre a face e o centroide das celulas utilizadas no stencil %
%
for i=1:face_num
    if i==1
        d1=0;
        d2=x(i)-faces(i);
    elseif i==face_num
        d1=0;
        d2=x(i-1)-faces(i);
    else
        d1=x(i)-faces(i);
        d2=x(i-1)-faces(i);
    end
    %
    % Valores dos Coeficientes do FDM %
    %
%     coefs(i,1)=1/(d1-d2);
%     coefs(i,2)=-1/(d1-d2);


    C=[ 1    1 ;
       d1   d2  ];
   
    e=[1;0];
    %
    const=inv(C)*e;

    coefs(i,1)=const(1);
    coefs(i,2)=const(2);


end
%
% Construção da Matriz %
%
% A=zeros(cell_num);
A=sparse([]);
b=sparse([]);
for i=1:cell_num
    if i==1
        A(i,i)=coefs(i+1,2)-coefs(i,2);
        A(i,i+1)=coefs(i+1,1);
%         b(i)=lap_phi(i)*cell_vol(i)-(0-coefs(i,1))*phi_West;
        b(i)=-(0-coefs(i,1))*phi_West;
    elseif i==cell_num
        A(i,i-1)=-coefs(i,2);
        A(i,i)=coefs(i+1,2)-coefs(i,1);
        b(i)=-coefs(i+1,1)*phi_East;
    else
        A(i,i-1)=-coefs(i,2);
        A(i,i)=coefs(i+1,2)-coefs(i,1);
        A(i,i+1)=coefs(i+1,1);
        b(i)=0;
    end
end
%
end
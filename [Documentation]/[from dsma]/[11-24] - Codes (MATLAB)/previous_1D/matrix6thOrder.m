function [A,b]=matrix6thOrder
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
    C=[ 1    1    1    1    1    1  ;
       d1   d2   d3   d4   d5   d6  ;
       d1^2 d2^2 d3^2 d4^2 d5^2 d6^2;
       d1^3 d2^3 d3^3 d4^3 d5^3 d6^3;
       d1^4 d2^4 d3^4 d4^4 d5^4 d6^4;
       d1^5 d2^5 d3^5 d4^5 d5^5 d6^5];
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
% Construção da Matriz %
%
% A=zeros(cell_num);
%A=sparse([]);
%b=sparse([]);
for i=1:cell_num
    if i==1
        A(i,i)=coefs(i+1,5)-coefs(i+0,2);
        A(i,i+1)=coefs(i+1,1)-coefs(i+0,3);
        A(i,i+2)=coefs(i+1,2)-coefs(i+0,4);
        A(i,i+3)=coefs(i+1,3)-coefs(i+0,5);
        A(i,i+4)=coefs(i+1,4)-coefs(i+0,6);
        b(i)=-(coefs(i+1,6)-coefs(i+0,1))*phi_West;
    elseif i==2
        A(i,i-1)=coefs(i+1,5)-coefs(i+0,5);
        A(i,i)=coefs(i+1,4)-coefs(i+0,1);
        A(i,i+1)=coefs(i+1,1)-coefs(i+0,2);
        A(i,i+2)=coefs(i+1,2)-coefs(i+0,3);
        A(i,i+3)=coefs(i+1,3)-coefs(i+0,4);
        b(i)=-(coefs(i+1,6)-coefs(i+0,6))*phi_West;
    elseif i==3
        A(i,i-2)=coefs(i+1,6)-coefs(i+0,5);
        A(i,i-1)=coefs(i+1,5)-coefs(i+0,4);
        A(i,i)=coefs(i+1,4)-coefs(i+0,1);
        A(i,i+1)=coefs(i+1,1)-coefs(i+0,2);
        A(i,i+2)=coefs(i+1,2)-coefs(i+0,3);
        A(i,i+3)=coefs(i+1,3);
        b(i)=+coefs(i+0,6)*phi_West;
    elseif i==cell_num-2
        A(i,i-3)=-coefs(i+0,6);
        A(i,i-2)=coefs(i+1,6)-coefs(i+0,5);
        A(i,i-1)=coefs(i+1,5)-coefs(i+0,4);
        A(i,i)=coefs(i+1,4)-coefs(i+0,1);
        A(i,i+1)=coefs(i+1,1)-coefs(i+0,2);
        A(i,i+2)=coefs(i+1,2)-coefs(i+0,3);
        b(i)=-coefs(i+1,3)*phi_East;
    elseif i==cell_num-1
        A(i,i-3)=coefs(i+1,6)-coefs(i+0,6);
        A(i,i-2)=coefs(i+1,5)-coefs(i+0,5);
        A(i,i-1)=coefs(i+1,4)-coefs(i+0,4);
        A(i,i)=coefs(i+1,3)-coefs(i+0,1);
        A(i,i+1)=coefs(i+1,1)-coefs(i+0,2);
        b(i)=-(coefs(i+1,2)-coefs(i+0,3))*phi_East;
    elseif i==cell_num
        A(i,i-4)=coefs(i+1,6)-coefs(i+0,6);
        A(i,i-3)=coefs(i+1,5)-coefs(i+0,5);
        A(i,i-2)=coefs(i+1,4)-coefs(i+0,4);
        A(i,i-1)=coefs(i+1,3)-coefs(i+0,3);
        A(i,i)=coefs(i+1,2)-coefs(i+0,1);
        b(i)=-(coefs(i+1,1)-coefs(i+0,2))*phi_East;
    else
        A(i,i-3)=-coefs(i+0,6);
        A(i,i-2)=coefs(i+1,6)-coefs(i+0,5);
        A(i,i-1)=coefs(i+1,5)-coefs(i+0,4);
        A(i,i)=coefs(i+1,4)-coefs(i+0,1);
        A(i,i+1)=coefs(i+1,1)-coefs(i+0,2);
        A(i,i+2)=coefs(i+1,2)-coefs(i+0,3);
        A(i,i+3)=coefs(i+1,3);
        b(i)=0;
    end
end
%
end
function [A,b]=matrix4thOrder
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
    % Valores dos Coeficientes do FDM %
    %
    coefs(i,1)=-(d2*d3 + d2*d4 + d3*d4)/(d1^2*d2 + d1^2*d3 + d1^2*d4 - d1^3 - d1*d2*d3 - d1*d2*d4 - d1*d3*d4 + d2*d3*d4);
    coefs(i,2)=-(d1*d3 + d1*d4 + d3*d4)/(d1*d2^2 + d2^2*d3 + d2^2*d4 - d2^3 - d1*d2*d3 - d1*d2*d4 + d1*d3*d4 - d2*d3*d4);
    coefs(i,3)=-(d1*d2 + d1*d4 + d2*d4)/(d1*d3^2 + d2*d3^2 + d3^2*d4 - d3^3 - d1*d2*d3 + d1*d2*d4 - d1*d3*d4 - d2*d3*d4);
    coefs(i,4)=-(d1*d2 + d1*d3 + d2*d3)/(d1*d4^2 + d2*d4^2 + d3*d4^2 - d4^3 + d1*d2*d3 - d1*d2*d4 - d1*d3*d4 - d2*d3*d4);
end
%
% Construção da Matriz %
%
% A=zeros(cell_num);
A=sparse([]);
b=sparse([]);
for i=1:cell_num
    if i==1
        A(i,i)=coefs(i+1,3)-coefs(i,2);
        A(i,i+1)=coefs(i+1,1)-coefs(i,3);
        A(i,i+2)=coefs(i+1,2)-coefs(i,4);
        b(i)=0+(coefs(i,1)-coefs(i+1,4))*phi_West;
    elseif i==2
        A(i,i-1)=coefs(i+1,4)-coefs(i,3);
        A(i,i)=coefs(i+1,3)-coefs(i,1);
        A(i,i+1)=coefs(i+1,1)-coefs(i,2);
        A(i,i+2)=coefs(i+1,2);
        b(i)=0+coefs(i,4)*phi_West;
    elseif i==cell_num-1
        A(i,i-2)=-coefs(i,4);
        A(i,i-1)=coefs(i+1,4)-coefs(i,3);
        A(i,i)=coefs(i+1,3)-coefs(i,1);
        A(i,i+1)=coefs(i+1,1)-coefs(i,2);
        b(i)=0-coefs(i+1,2)*phi_East;
    elseif i==cell_num
        A(i,i-2)=coefs(i+1,4)-coefs(i,4);
        A(i,i-1)=coefs(i+1,3)-coefs(i,3);
        A(i,i)=coefs(i+1,2)-coefs(i,1);
        b(i)=0+(coefs(i,2)-coefs(i+1,1))*phi_East;
    else
        A(i,i-2)=-coefs(i,4);
        A(i,i-1)=coefs(i+1,4)-coefs(i,3);
        A(i,i)=coefs(i+1,3)-coefs(i,1);
        A(i,i+1)=coefs(i+1,1)-coefs(i,2);
        A(i,i+2)=coefs(i+1,2);
        b(i)=0;
    end
end
%    
end
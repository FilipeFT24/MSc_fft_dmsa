function f=solutions(solution,coord,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   07 de Junho de 2016                                             %
%                                                                                                   %
% Função que determina as soluções analiticas                                                       %
%                                                                                                   %
% sin - Função sinuisoidal;                                                                         %
% exp - Função exponencial;                                                                         %
%                                                                                                   %
% anal  - Solução Analitica;                                                                        %
% lap   - Laplaciano (2ª Derivada);                                                                 %
% flux  - Fluxo (1ª Derivada);                                                                      %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
global L Lref cell_num face_num;
global x faces cell_vol;
global u_convec gamma_diff;
%
if strcmp(solution,'sin')==1
    if strcmp(type,'anal')==1
        f=sin(3*pi*(coord));
    end
    if strcmp(type,'lap')==1
        f=-(3*pi)^2*sin(3*pi*(coord));
    end
    if strcmp(type,'flux')==1
        f=3*pi*cos(3*pi*(coord));
    end
end
if strcmp(solution,'exp')==1
    %
    s=sqrt(0.0175);
    u=0.5;
    %
    if strcmp(type,'anal')==1
        f=exp(-(coord-u)^2/s^2);
    end
    if strcmp(type,'lap')==1
%         f=-2/s^2*(1-2*(coord-u)^2/s^2)*exp(-(coord-u)^2/s^2);
        f=((4*coord^2-8*u*coord+4*u^2-2*s^2)*exp((-(coord-u)^2/s^2)))/s^4;
    end
    if strcmp(type,'flux')==1
        f=-2*((coord-u)/s^2)*exp(-(coord-u)^2/s^2);
    end
end
if strcmp(solution,'zero')==1
    
    if strcmp(type,'anal')==1
%         f=1+(1-exp(25*coord))/(7.2*10^10);
%           f=(2.7183-exp(coord))/1.7183;
            f=(1-exp(-u_convec/gamma_diff*coord)/(exp(-1*u_convec/gamma_diff)));
    end
    if strcmp(type,'lap')==1
%         f=-2/s^2*(1-2*(coord-u)^2/s^2)*exp(-(coord-u)^2/s^2);
        f=0;
    end
    if strcmp(type,'flux')==1
        f=0;
    end
end
if strcmp(solution,'x')==1
    
    if strcmp(type,'anal')==1
            f= ((u_convec^2 + 4*gamma_diff*coord)^(3/2)/(6*gamma_diff) - u_convec*coord)/(2*gamma_diff);
    end
    if strcmp(type,'lap')==1
        f= 1/sqrt(u_convec^2 + 4*gamma_diff*coord);
    end
    if strcmp(type,'flux')==1
        f= (sqrt(u_convec^2 + 4*gamma_diff*coord) - u_convec)/(2*gamma_diff);%coord/u_convec;%
    end
end
%%
end
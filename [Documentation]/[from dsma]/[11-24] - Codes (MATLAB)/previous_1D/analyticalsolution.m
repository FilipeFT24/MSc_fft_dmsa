function [phi,lap_phi,phi_0,phi_L,phi_flux_0,phi_flux_L]=analyticalsolution(solution,metodo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   07 de Junho de 2016                                             %
%                                                                                                   %
% Função que determina as soluções analiticas                                                       %
%                                                                                                   %
% phi           - Valor da Propriedade no centroide da célula;                                      %
% lap_phi       - Laplaciano da Propriedade na célula, utiliza pontos de Gauss;                     %
% phi_0         - Valor da Propriedade na fronteira esquerda;                                       %
% phi_L         - Valor da Propriedade na fronteira direita;                                        %
% phi_flux_0    - Valor do Fluxo na fronteira esquerda;                                             %
% phi_flux_L    - Valor do Fluxo na fronteira direira.                                              %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global u_convec gamma_diff;
global L Lref cell_num face_num;
global x faces cell_vol;
%
if metodo=='2ndOrder'
    for i=1:cell_num
        phi(i)=solutions(solution,x(i),'anal');
        % Calula o Laplaciano de Phi com base nos Pontos de Gauss %
        lap_phi(i)=gamma_diff*solutions(solution,x(i),'lap')+u_convec*solutions(solution,x(i),'flux');
    end
elseif metodo=='4thOrder'
    for i=1:cell_num
        x1=x(i);
        phi(i)=solutions(solution,x1,'anal');
        %Calula o Laplaciano de Phi com base nos Pontos de Gauss %
        x1=x(i)-cell_vol(i)/2*sqrt(1/3);
        x2=x(i)+cell_vol(i)/2*sqrt(1/3);
        lap1=gamma_diff*solutions(solution,x1,'lap')+u_convec*solutions(solution,x1,'flux');
        lap2=gamma_diff*solutions(solution,x2,'lap')+u_convec*solutions(solution,x2,'flux');
        lap_phi(i)=(lap1+lap2)/2;
    end
    s=1;
elseif metodo=='6thOrder'
    for i=1:cell_num
        phi(i)=solutions(solution,x(i),'anal');
       % Calula o Laplaciano de Phi com base nos Pontos de Gauss %
        x1=x(i)-cell_vol(i)/2*sqrt(3/5);
        x2=x(i)+cell_vol(i)/2*sqrt(3/5);
        lap1=gamma_diff*solutions(solution,x1,'lap')+u_convec*solutions(solution,x1,'flux');
        lap2=gamma_diff*solutions(solution,x2,'lap')+u_convec*solutions(solution,x2,'flux');
        lap3=gamma_diff*solutions(solution,x(i),'lap')+u_convec*solutions(solution,x(i),'flux');
        lap_phi(i)=(lap1+lap2)*5/18+lap3*4/9;
    end 
elseif metodo=='8thOrder'
    for i=1:cell_num
        phi(i)=solutions(solution,x(i),'anal');
        %Calula o Laplaciano de Phi com base nos Pontos de Gauss %
        x1=x(i)-cell_vol(i)/2*1/35*sqrt(525-70*sqrt(30));
        x2=x(i)-cell_vol(i)/2*1/35*sqrt(525+70*sqrt(30));
        x3=x(i)+cell_vol(i)/2*1/35*sqrt(525-70*sqrt(30));
        x4=x(i)+cell_vol(i)/2*1/35*sqrt(525+70*sqrt(30));
        lap1=gamma_diff*solutions(solution,x1,'lap')+u_convec*solutions(solution,x1,'flux');
        lap2=gamma_diff*solutions(solution,x2,'lap')+u_convec*solutions(solution,x2,'flux');
        lap3=gamma_diff*solutions(solution,x3,'lap')+u_convec*solutions(solution,x3,'flux');
        lap4=gamma_diff*solutions(solution,x4,'lap')+u_convec*solutions(solution,x4,'flux');
        lap_phi(i)=(lap1+lap3)*(18+sqrt(30))/(2*36)+(lap2+lap4)*(18-sqrt(30))/(2*36);
    end
end
% Valores na Fronteira %
phi_0=solutions(solution,0,'anal');
phi_L=solutions(solution,L,'anal');
phi_flux_0=solutions(solution,0,'flux');
phi_flux_L=solutions(solution,L,'flux');
%
end
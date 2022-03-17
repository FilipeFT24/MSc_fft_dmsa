function [phi,lap_phi]=analyticalsourceterm(solution,metodo)
%
global L Lref cell_num face_num;
global x faces cell_vol;
%
if solution=='sin'
    phi=sin(3*pi*x);
    lap_phi=-(3*pi)^2*sin(3*pi*x);
elseif solution=='exp'
    %
    s=sqrt(0.0175);
    u=0.5;
    %
    phi=exp(-(x-u)^2/s^2);
    lap_phi=-2/s^2*(1-2*(x-u)^2/s^2)*exp(-(x-u)^2/s^2);
else
    error('Solução Analitica Não Implementada');
end 
end
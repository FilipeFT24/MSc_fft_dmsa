function f=lap_solutions(solution,coord)
%
global L Lref cell_num face_num;
global x faces cell_vol;
%
if solution=='sin'
    f=-(3*pi)^2*sin(3*pi*coord);    
end
if solution=='exp'
    %
    s=sqrt(0.0175);
    u=0.5;
    %
    f=-2/s^2*(1-2*(coord-u)^2/s^2)*exp(-(coord-u)^2/s^2);
end
end

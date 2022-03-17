function f=PolyReconstruction6thOrder(C,x,y,xc,yc,type,face)

global L Lref cell_side vert_side cell_num face_num vert_num w;
global verts cells cell_verts cell_faces cell_vol faces face_area face_vert cell_norm;
global face_bound cell_bound face_cells vert_cells;
global phi lap_phi  phi_faces flux_phi_faces face_w_gauss extended_stencil exact_coeffs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                  14 de Outubro de 2016                                            %
%                                  14 de Outubro de 2016                                            %
%                                                                                                   %
% Fun��o calcula a contirbui��o de cada celula para a reconstru��o do polinomio na face que se esta %
% a trabalhar, esta fun��o tem como objectivo facilitar a constru��o da matriz A.                   %
% Aten��o que esta fun��o n�o determina o polinomio reconstruido na face, serve apenas para quando  %
% se esta a determinar as entradas na matriz A.                                                     %
% Caso queira saber o polinomio da face tenho de usar outra fun��o                                  %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  % Normal Exterior da Face %
    vert=face_vert(face,:);
    nx=verts(vert(1),1)-verts(vert(2),1);
    ny=verts(vert(1),2)-verts(vert(2),2);



if strcmp(type,'poly')==1
    if extended_stencil && nx==0
        f=C(1)+C(2)*(x-xc)+C(3)*(y-yc) ...
        +C(4)*(x-xc)^2+C(5)*(y-yc)^2+C(6)*(x-xc)*(y-yc) ...
        +C(7)*(x-xc)^3+C(8)*(y-yc)^3+C(9)*(x-xc)^2*(y-yc)+C(10)*(x-xc)*(y-yc)^2 ...
        +C(11)*(x-xc)^4+C(12)*(y-yc)^4+C(13)*(x-xc)^3*(y-yc)+C(14)*(x-xc)^2*(y-yc)^2+C(15)*(x-xc)*(y-yc)^3 ...
        +C(16)*(x-xc)^5+C(17)*(y-yc)^5+C(18)*(x-xc)^4*(y-yc)+C(19)*(x-xc)^3*(y-yc)^2+C(20)*(x-xc)^2*(y-yc)^3+C(21)*(x-xc)*(y-yc)^4 ...
        +C(22)*(y-yc)^6+C(23)*(x-xc)^5*(y-yc)+C(24)*(x-xc)^4*(y-yc)^2+C(25)*(x-xc)^3*(y-yc)^3+C(26)*(x-xc)^2*(y-yc)^4+C(27)*(x-xc)^1*(y-yc)^5;
    elseif extended_stencil && ny==0
        f=C(1)+C(2)*(x-xc)+C(3)*(y-yc) ...
        +C(4)*(x-xc)^2+C(5)*(y-yc)^2+C(6)*(x-xc)*(y-yc) ...
        +C(7)*(x-xc)^3+C(8)*(y-yc)^3+C(9)*(x-xc)^2*(y-yc)+C(10)*(x-xc)*(y-yc)^2 ...
        +C(11)*(x-xc)^4+C(12)*(y-yc)^4+C(13)*(x-xc)^3*(y-yc)+C(14)*(x-xc)^2*(y-yc)^2+C(15)*(x-xc)*(y-yc)^3 ...
        +C(16)*(x-xc)^5+C(17)*(y-yc)^5+C(18)*(x-xc)^4*(y-yc)+C(19)*(x-xc)^3*(y-yc)^2+C(20)*(x-xc)^2*(y-yc)^3+C(21)*(x-xc)*(y-yc)^4 ...
        +C(22)*(x-xc)^6+C(23)*(x-xc)^5*(y-yc)+C(24)*(x-xc)^4*(y-yc)^2+C(25)*(x-xc)^3*(y-yc)^3+C(26)*(x-xc)^2*(y-yc)^4+C(27)*(x-xc)^1*(y-yc)^5;
   
    elseif exact_coeffs
        f=C(1)+C(2)*(x-xc)+C(3)*(y-yc) ...
        +C(4)*(x-xc)^2+C(5)*(y-yc)^2+C(6)*(x-xc)*(y-yc) ...
        +C(7)*(x-xc)^3+C(8)*(y-yc)^3+C(9)*(x-xc)^2*(y-yc)+C(10)*(x-xc)*(y-yc)^2 ...
        +C(11)*(x-xc)^4+C(12)*(y-yc)^4+C(13)*(x-xc)^3*(y-yc)+C(14)*(x-xc)^2*(y-yc)^2+C(15)*(x-xc)*(y-yc)^3 ...
        +C(16)*(x-xc)^5+C(17)*(y-yc)^5+C(18)*(x-xc)^4*(y-yc)+C(19)*(x-xc)^3*(y-yc)^2+C(20)*(x-xc)^2*(y-yc)^3+C(21)*(x-xc)*(y-yc)^4 ...
        +C(22)*(x-xc)^5*(y-yc)+C(23)*(x-xc)^4*(y-yc)^2+C(24)*(x-xc)^3*(y-yc)^3+C(25)*(x-xc)^2*(y-yc)^4+C(26)*(x-xc)^1*(y-yc)^5 ...
        +C(27)*(x-xc)^5*(y-yc)^2+C(28)*(x-xc)^4*(y-yc)^3+C(29)*(x-xc)^3*(y-yc)^4+C(30)*(x-xc)^2*(y-yc)^5 ... 
        +C(31)*(x-xc)^5*(y-yc)^3+C(32)*(x-xc)^4*(y-yc)^4+C(33)*(x-xc)^3*(y-yc)^5 ... 
        +C(34)*(x-xc)^5*(y-yc)^4+C(35)*(x-xc)^4*(y-yc)^5 ...
        +C(36)*(x-xc)^5*(y-yc)^5;
    
    if nx==0  
     f=f+C(37)*(y-yc)^6+C(38)*(y-yc)^6*(x-xc)+C(39)*(x-xc)^2*(y-yc)^6+C(40)*(x-xc)^3*(y-yc)^6+C(41)*(x-xc)^4*(y-yc)^6+C(42)*(x-xc)^5*(y-yc)^6;   
    end
    
    if ny==0  
     f=f+C(37)*(x-xc)^6+C(38)*(y-yc)^1*(x-xc)^6+C(39)*(x-xc)^6*(y-yc)^2+C(40)*(x-xc)^6*(y-yc)^3+C(41)*(x-xc)^6*(y-yc)^4+C(42)*(x-xc)^6*(y-yc)^5;   
    end
        
      
    
    else
    f=C(1)+C(2)*(x-xc)+C(3)*(y-yc) ...
        +C(4)*(x-xc)^2+C(5)*(y-yc)^2+C(6)*(x-xc)*(y-yc) ...
        +C(7)*(x-xc)^3+C(8)*(y-yc)^3+C(9)*(x-xc)^2*(y-yc)+C(10)*(x-xc)*(y-yc)^2 ...
        +C(11)*(x-xc)^4+C(12)*(y-yc)^4+C(13)*(x-xc)^3*(y-yc)+C(14)*(x-xc)^2*(y-yc)^2+C(15)*(x-xc)*(y-yc)^3 ...
        +C(16)*(x-xc)^5+C(17)*(y-yc)^5+C(18)*(x-xc)^4*(y-yc)+C(19)*(x-xc)^3*(y-yc)^2+C(20)*(x-xc)^2*(y-yc)^3+C(21)*(x-xc)*(y-yc)^4;
    end
elseif strcmp(type,'xflux')==1
     if extended_stencil && nx==0
        f=C(2) ...
        +2*C(4)*(x-xc)+C(6)*(y-yc) ...
        +3*C(7)*(x-xc)^2+2*C(9)*(x-xc)*(y-yc)+C(10)*(y-yc)^2 ...
        +4*C(11)*(x-xc)^3+3*C(13)*(x-xc)^2*(y-yc)+2*C(14)*(x-xc)*(y-yc)^2+C(15)*(y-yc)^3 ...
        +5*C(16)*(x-xc)^4+4*C(18)*(x-xc)^3*(y-yc)+3*C(19)*(x-xc)^2*(y-yc)^2+2*C(20)*(x-xc)*(y-yc)^3+C(21)*(y-yc)^4 ...
        +5*C(23)*(x-xc)^4*(y-yc)+4*C(24)*(x-xc)^3*(y-yc)^2+3*C(25)*(x-xc)^2*(y-yc)^3+2*C(26)*(x-xc)^1*(y-yc)^4+C(27)*(y-yc)^5;
    elseif extended_stencil && ny==0
        f=C(2) ...
        +2*C(4)*(x-xc)+C(6)*(y-yc) ...
        +3*C(7)*(x-xc)^2+2*C(9)*(x-xc)*(y-yc)+C(10)*(y-yc)^2 ...
        +4*C(11)*(x-xc)^3+3*C(13)*(x-xc)^2*(y-yc)+2*C(14)*(x-xc)*(y-yc)^2+C(15)*(y-yc)^3 ...
        +5*C(16)*(x-xc)^4+4*C(18)*(x-xc)^3*(y-yc)+3*C(19)*(x-xc)^2*(y-yc)^2+2*C(20)*(x-xc)*(y-yc)^3+C(21)*(y-yc)^4 ...
        +6*C(22)*(x-xc)^5+5*C(23)*(x-xc)^4*(y-yc)+4*C(24)*(x-xc)^3*(y-yc)^2+3*C(25)*(x-xc)^2*(y-yc)^3+2*C(26)*(x-xc)^1*(y-yc)^4+C(27)*(y-yc)^5;
    
     elseif exact_coeffs
        f=C(2) ...
        +2*C(4)*(x-xc)+C(6)*(y-yc) ...
        +3*C(7)*(x-xc)^2+2*C(9)*(x-xc)*(y-yc)+C(10)*(y-yc)^2 ...
        +4*C(11)*(x-xc)^3+3*C(13)*(x-xc)^2*(y-yc)+2*C(14)*(x-xc)*(y-yc)^2+C(15)*(y-yc)^3 ...
        +5*C(16)*(x-xc)^4+4*C(18)*(x-xc)^3*(y-yc)+3*C(19)*(x-xc)^2*(y-yc)^2+2*C(20)*(x-xc)*(y-yc)^3+C(21)*(y-yc)^4 ...
        +5*C(22)*(x-xc)^4*(y-yc)+4*C(23)*(x-xc)^3*(y-yc)^2+3*C(24)*(x-xc)^2*(y-yc)^3+2*C(25)*(x-xc)^1*(y-yc)^4+C(26)*(y-yc)^5 ...
        +5*C(27)*(x-xc)^4*(y-yc)^2+4*C(28)*(x-xc)^3*(y-yc)^3+3*C(29)*(x-xc)^2*(y-yc)^4+2*C(30)*(x-xc)^1*(y-yc)^5 ... 
        +5*C(31)*(x-xc)^4*(y-yc)^3+4*C(32)*(x-xc)^3*(y-yc)^4+3*C(33)*(x-xc)^2*(y-yc)^5 ... 
        +5*C(34)*(x-xc)^4*(y-yc)^4+4*C(35)*(x-xc)^3*(y-yc)^5 ...
        +5*C(36)*(x-xc)^4*(y-yc)^5;
    
    if nx==0  
     f=f+C(38)*(y-yc)^6+2*C(39)*(x-xc)^1*(y-yc)^6+3*C(40)*(x-xc)^2*(y-yc)^6+4*C(41)*(x-xc)^3*(y-yc)^6+5*C(42)*(x-xc)^4*(y-yc)^6;   
    end
    
    if ny==0  
     f=f+6*C(37)*(x-xc)^5+6*C(38)*(y-yc)^1*(x-xc)^5+6*C(39)*(x-xc)^5*(y-yc)^2+6*C(40)*(x-xc)^5*(y-yc)^3+6*C(41)*(x-xc)^5*(y-yc)^4+6*C(42)*(x-xc)^5*(y-yc)^5;   
    end
     
     
     
     
    else
    f=C(2) ...
        +2*C(4)*(x-xc)+C(6)*(y-yc) ...
        +3*C(7)*(x-xc)^2+2*C(9)*(x-xc)*(y-yc)+C(10)*(y-yc)^2 ...
        +4*C(11)*(x-xc)^3+3*C(13)*(x-xc)^2*(y-yc)+2*C(14)*(x-xc)*(y-yc)^2+C(15)*(y-yc)^3 ...
        +5*C(16)*(x-xc)^4+4*C(18)*(x-xc)^3*(y-yc)+3*C(19)*(x-xc)^2*(y-yc)^2+2*C(20)*(x-xc)*(y-yc)^3+C(21)*(y-yc)^4;
     end
elseif strcmp(type,'yflux')==1
     if extended_stencil && nx==0
        f=C(3) ...
        +2*C(5)*(y-yc)+C(6)*(x-xc) ...
        +3*C(8)*(y-yc)^2+C(9)*(x-xc)^2+2*C(10)*(x-xc)*(y-yc) ...
        +4*C(12)*(y-yc)^3+C(13)*(x-xc)^3+2*C(14)*(x-xc)^2*(y-yc)+3*C(15)*(x-xc)*(y-yc)^2 ...
        +5*C(17)*(y-yc)^4+C(18)*(x-xc)^4+2*C(19)*(x-xc)^3*(y-yc)+3*C(20)*(x-xc)^2*(y-yc)^2+4*C(21)*(x-xc)*(y-yc)^3 ...
        +6*C(22)*(y-yc)^5+C(23)*(x-xc)^5+2*C(24)*(x-xc)^4*(y-yc)^1+3*C(25)*(x-xc)^3*(y-yc)^2+4*C(26)*(x-xc)^2*(y-yc)^3+5*C(27)*(x-xc)^1*(y-yc)^4;
    elseif extended_stencil && ny==0
            f=C(3) ...
        +2*C(5)*(y-yc)+C(6)*(x-xc) ...
        +3*C(8)*(y-yc)^2+C(9)*(x-xc)^2+2*C(10)*(x-xc)*(y-yc) ...
        +4*C(12)*(y-yc)^3+C(13)*(x-xc)^3+2*C(14)*(x-xc)^2*(y-yc)+3*C(15)*(x-xc)*(y-yc)^2 ...
        +5*C(17)*(y-yc)^4+C(18)*(x-xc)^4+2*C(19)*(x-xc)^3*(y-yc)+3*C(20)*(x-xc)^2*(y-yc)^2+4*C(21)*(x-xc)*(y-yc)^3 ...
        +C(23)*(x-xc)^5+2*C(24)*(x-xc)^4*(y-yc)^1+3*C(25)*(x-xc)^3*(y-yc)^2+4*C(26)*(x-xc)^2*(y-yc)^3+5*C(27)*(x-xc)^1*(y-yc)^4;
   
     elseif exact_coeffs
        f=C(3) ...
        +2*C(5)*(y-yc)+C(6)*(x-xc) ...
        +3*C(8)*(y-yc)^2+C(9)*(x-xc)^2+2*C(10)*(x-xc)*(y-yc) ...
        +4*C(12)*(y-yc)^3+C(13)*(x-xc)^3+2*C(14)*(x-xc)^2*(y-yc)+3*C(15)*(x-xc)*(y-yc)^2 ...
        +5*C(17)*(y-yc)^4+C(18)*(x-xc)^4+2*C(19)*(x-xc)^3*(y-yc)+3*C(20)*(x-xc)^2*(y-yc)^2+4*C(21)*(x-xc)*(y-yc)^3 ...
        +C(22)*(x-xc)^5+2*C(23)*(x-xc)^4*(y-yc)^1+3*C(24)*(x-xc)^2*(y-yc)^3+4*C(25)*(x-xc)^2*(y-yc)^3+5*C(26)*(x-xc)^1*(y-yc)^4 ...
        +2*C(27)*(x-xc)^5*(y-yc)^1+3*C(28)*(x-xc)^4*(y-yc)^2+4*C(29)*(x-xc)^3*(y-yc)^3+5*C(30)*(x-xc)^2*(y-yc)^4 ... 
        +3*C(31)*(x-xc)^5*(y-yc)^2+4*C(32)*(x-xc)^4*(y-yc)^3+5*C(33)*(x-xc)^3*(y-yc)^4 ... 
        +4*C(34)*(x-xc)^5*(y-yc)^3+5*C(35)*(x-xc)^4*(y-yc)^4 ...
        +5*C(36)*(x-xc)^5*(y-yc)^4;
    
    if nx==0  
     f=f+6*C(37)*(y-yc)^5+6*C(38)*(y-yc)^5*(x-xc)+6*C(39)*(x-xc)^2*(y-yc)^5+6*C(40)*(x-xc)^3*(y-yc)^5+6*C(41)*(x-xc)^4*(y-yc)^5+6*C(42)*(x-xc)^5*(y-yc)^5;   
    end
    
    if ny==0  
     f=f+C(38)*(x-xc)^6+2*C(39)*(x-xc)^6*(y-yc)^1+3*C(40)*(x-xc)^6*(y-yc)^2+4*C(41)*(x-xc)^6*(y-yc)^3+5*C(42)*(x-xc)^6*(y-yc)^4;   
    end
     
     
     
     
     
    
     else
    f=C(3) ...
        +2*C(5)*(y-yc)+C(6)*(x-xc) ...
        +3*C(8)*(y-yc)^2+C(9)*(x-xc)^2+2*C(10)*(x-xc)*(y-yc) ...
        +4*C(12)*(y-yc)^3+C(13)*(x-xc)^3+2*C(14)*(x-xc)^2*(y-yc)+3*C(15)*(x-xc)*(y-yc)^2 ...
        +5*C(17)*(y-yc)^4+C(18)*(x-xc)^4+2*C(19)*(x-xc)^3*(y-yc)+3*C(20)*(x-xc)^2*(y-yc)^2+4*C(21)*(x-xc)*(y-yc)^3;
     end
else
    error('\n\n\tERRO: Fluxo Desconhecido\n\n');
end
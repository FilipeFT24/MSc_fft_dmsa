function f=SolutionDiffusion(solution,x,y,type)
global u_convec_x u_convec_y gamma_diff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   28 de Junho de 2016                                             %
%                                  12 de Outubro de 2016                                            %
%                                                                                                   %
% Fun��o que determina as solu��es analiticas                                                       %
%                                                                                                   %
% solution - Fun��o Analitica que se pretende calcular                                              %
% x        - Coordenada X                                                                           %
% y        - Coordenada Y                                                                           %
% type     - Escolhe entre a solu��o analitica 'anal', laplaciano 'lap' e as derivadas em x e y     %
%            'xflux' ou 'yflux'                                                                     %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmp(solution,'sin')==1
    %
    % Solu��o Analitica %
    if strcmp(type,'anal')==1
        f=sin(3*pi*x)*sin(3*pi*y);
    end
    %
    % Laplaciano %
    if strcmp(type,'lap')==1
        f=-18*pi^2*sin(3*pi*x)*sin(3*pi*y);
    end
    %
    % Derivada X %
    if strcmp(type,'xflux')==1
        f=3*pi*cos(3*pi*x)*sin(3*pi*y);
    end
    %
    % Derivada Y %
    if strcmp(type,'yflux')==1
        f=3*pi*sin(3*pi*x)*cos(3*pi*y);
    end

elseif strcmp(solution,'sin2')==1
    %
    % Solu��o Analitica %
    if strcmp(type,'anal')==1
        f=50000*cos(3*pi*x)*sin(2*pi*y);
    end
    %
    % Laplaciano %
    if strcmp(type,'lap')==1
        f=-9*50000*pi^2*cos(3*pi*x)*sin(2*pi*y)-50000*4*pi^2*cos(3*pi*x)*sin(2*pi*y);
    end
    %
    % Derivada X %
    if strcmp(type,'xflux')==1
        f=-3*50000*pi*sin(3*pi*x)*sin(2*pi*y);
    end
    %
    % Derivada Y %
    if strcmp(type,'yflux')==1
        f=2*50000*pi*cos(3*pi*x)*cos(2*pi*y);
    end
    
    
elseif strcmp(solution,'exp')==1
    %
    % Constantes %
    s=sqrt(0.0175);
    u=0.5;
    %
    % Solu��o Analitica %
    if strcmp(type,'anal')==1
        f=exp(-((x-u)^2+(y-u)^2)/s^2);
    end
    %
    % Laplaciano %
    if strcmp(type,'lap')==1
        f=(-2/s^2+4*(x-u)^2/s^4+4*(y-u)^2/s^4-2/s^2)*exp(-((x-u)^2+(y-u)^2)/s^2);
        %-4/s^2*(1-((x-u)^2+(y-u)^2)/s^2)*exp(((x-u)^2+(y-u)^2)/s^2);
    end
    %
    % Derivada X %
    if strcmp(type,'xflux')==1
        f=-2*(x-u)/s^2*exp(-((x-u)^2+(y-u)^2)/s^2);
    end
    %
    % Derivada Y %
    if strcmp(type,'yflux')==1
        f=-2*(y-u)/s^2*exp(-((x-u)^2+(y-u)^2)/s^2);
    end
    
    
%% Em teste    
elseif strcmp(solution,'paper')==1
    %
    % Constantes %
    
    C = 65;%11236%65
    u = u_convec_x;
    v = u_convec_y;
    k = -gamma_diff; %difusivo
    
    % fun�oes auxiliares %
    
    alpha =@(x)1/u*(x-(exp(u*x)-1)/(exp(u)-1));
    beta =@(y)1/v*(y-(exp(v*y)-1)/(exp(v)-1));
    
 
    % Solu��o Analitica %
    if strcmp(type,'anal')==1
        f=C*alpha(x)*beta(y);
    end
    %
    % Laplaciano %
    if strcmp(type,'lap')==1
        %f=-C*(alpha(x)+beta(y));
        f=((1-u*exp(u*x)/(exp(u)-1))*C*beta(y)+(1-v*exp(v*y)/(exp(v)-1))*C*alpha(x)+(u*exp(u*x)/(exp(u)-1))*C*beta(y)*k+(v*exp(v*y)/(exp(v)-1))*C*alpha(x)*k);
    end
    %
    % Derivada X %
    if strcmp(type,'xflux')==1
        f=0;
    end
    %
    % Derivada Y %
    if strcmp(type,'yflux')==1
        f=0;
    end
%% FIm teste



elseif strcmp(solution,'2nd')==1
    %
    % Solu��o Analitica %
    if strcmp(type,'anal')==1
        f=1/2*(x+y);
    end
    %
    % Laplaciano %
    if strcmp(type,'lap')==1
        f=0;
    end
    %
    % Derivada X %
    if strcmp(type,'xflux')==1
        f=1/2;
    end
    %
    % Derivada Y %
    if strcmp(type,'yflux')==1
        f=1/2;
    end
elseif strcmp(solution,'4th')==1
    %
    % Solu��o Analitica %
    if strcmp(type,'anal')==1
        f=(y-1)^3+(x-1)^3;
    end
    %
    % Laplaciano %
    if strcmp(type,'lap')==1
        f=6*(y-1)+6*(x-1);
    end
    %
    % Derivada X %
    if strcmp(type,'xflux')==1
        f=3*(x-1)^2;
    end
    %
    % Derivada Y %
    if strcmp(type,'yflux')==1
        f=3*(y-1)^2;
    end
elseif strcmp(solution,'6th')==1
    %
    % Solu��o Analitica %
    if strcmp(type,'anal')==1
        f=x^5;
    end
    %
    % Laplaciano %
    if strcmp(type,'lap')==1
        f=20*x^3;
    end
    %
    % Derivada X %
    if strcmp(type,'xflux')==1
        f=5*x^4;
    end
    %
    % Derivada Y %
    if strcmp(type,'yflux')==1
        f=0;
    end
elseif strcmp(solution,'8th')==1
    %
    % Solu��o Analitica %
    if strcmp(type,'anal')==1
        f=y^3;
    end
    %
    % Laplaciano %
    if strcmp(type,'lap')==1
        f=6*y;
    end
    %
    % Derivada X %
    if strcmp(type,'xflux')==1
        f=0;
    end
    %
    % Derivada Y %
    if strcmp(type,'yflux')==1
        f=3*y^2;
    end
    
 elseif strcmp(solution,'9th')==1
    %
    % Solu��o Analitica %
    if strcmp(type,'anal')==1
        f=(x^7)*(y^8) + (x^8)*(y^7);
    end
    %
    % Laplaciano %
    if strcmp(type,'lap')==1
        f=42*(x^5)*(y^8)+56*(x^6)*(y^7)+56*(x^7)*(y^6)+42*(x^8)*(y^5);
    end
    %
    % Derivada X %
    if strcmp(type,'xflux')==1
        f=7*(x^6)*(y^8)+8*(x^7)*(y^7);
    end
    %
    % Derivada Y %
    if strcmp(type,'yflux')==1
        f=8*(x^7)*(y^7)+7*(x^8)*(y^6);
    end
    
 elseif strcmp(solution,'10th')==1
    %
    % Solu��o Analitica %
    if strcmp(type,'anal')==1
        f=(x^7)*(y^7);
    end
    %
    % Laplaciano %
    if strcmp(type,'lap')==1
        f=42*(x^5)*(y^7)+42*(x^7)*(y^5);
    end
    %
    % Derivada X %
    if strcmp(type,'xflux')==1
        f=7*(x^6)*(y^7);
    end
    %
    % Derivada Y %
    if strcmp(type,'yflux')==1
        f=7*(x^7)*(y^6);
    end    
 elseif strcmp(solution,'paper2')==1
    %
    % Solu��o Analitica %
    if strcmp(type,'anal')==1
        f=4*y*(1-y)/(x+1);
    end
    %
    % Laplaciano %
    if strcmp(type,'lap')==1
        f=-(8*(y-1)*y/(x+1)^3)-8/(x+1);
    end
    %
    % Derivada X %
    if strcmp(type,'xflux')==1
        f=4*(y-1)*y/(x+1)^2;
    end
    %
    % Derivada Y %
    if strcmp(type,'yflux')==1
        f=(4-8*y)/(x+1);
    end
elseif strcmp(solution,'11th')==1
    %
    % Solu��o Analitica %
    if strcmp(type,'anal')==1
        f=exp(x);
    end
    %
    % Laplaciano %
    if strcmp(type,'lap')==1
        f=exp(x);
    end
    %
    % Derivada X %
    if strcmp(type,'xflux')==1
        f=exp(x);
    end
    %
    % Derivada Y %
    if strcmp(type,'yflux')==1
        f=0;
    end       
    
elseif strcmp(solution,'poly3')==1
    %
    % Solu��o Analitica %
    if strcmp(type,'anal')==1
        f=1+1*x+1*y+1*x^2+1*y^2+1*x*y+1*x^3+1*y^3;%+1*x^2*y+1*y^2*x;
    end
    %
    % Laplaciano %
    if strcmp(type,'lap')==1
        f=2+6*x+2+6*y;%+2*x+2*y;
    end
    %
    % Derivada X %
    if strcmp(type,'xflux')==1
        f= 1+2*x+y+3*x^2;%+2*x*y+y^2;
    end
    %
    % Derivada Y %
    if strcmp(type,'yflux')==1
        f= 1+2*y+x+3*y^2;%+x^2+2*x*y;
    end       
end
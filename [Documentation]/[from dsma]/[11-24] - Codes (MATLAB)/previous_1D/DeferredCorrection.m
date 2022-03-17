function [phi_num,iter,res]=DeferredCorrection(metodo1,metodo2)
%
%% Declaração das Variaveis Globais %%
%
global u_convec gamma_diff numero_condicao fl0;
global L Lref cell_num face_num;
global x faces cell_vol;
global phi lap_phi phi_West phi_East;
%
%% Construção das Matrizes do Metodo %%
%
if metodo1=='2ndOrder'
    [Ad,bd]=matrix2ndOrder;  %difusão
    [Ac,bc]=matrix2ndOrder_convection; %convecçao
elseif metodo1=='4thOrder'
    [Ad,bd]=matrix4thOrder;
    [Ac,bc]=matrix4thOrder_convection;
elseif metodo1=='6thOrder'
    [Ad,bd]=matrix6thOrder;
    [Ac,bc]=matrix6thOrder_convection;
elseif metodo1=='8thOrder'
    [Ad,bd]=matrix8thOrder;
    [Ac,bc]=matrix8thOrder_convection;
end

[b]=matrix_b;
b1 =b+gamma_diff*bd+u_convec*bc;

% if metodo2=='2ndOrder'
%     [A2,b2]=matrix2ndOrder;
% elseif metodo2=='4thOrder'
%     [A2,b2]=matrix4thOrder;
% elseif metodo2=='6thOrder'
%     [A2,b2]=matrix6thOrder;
% elseif metodo2=='8thOrder'
%     [A2,b2]=matrix8thOrder;
% end
%
%% Método da Correção Diferida %%
%


sol_ini = 0*ones(size(phi,2),1);
I = sparse(eye(size(Ac,1)));

% Caso se utilze correcçao diferida, este facto pode ajudar a estabilizar
% Caso não se utilize, deve estar com o valor zero

factor = 0;%diagonal_dominance(gamma_diff*Ad+u_convec*Ac);

%inicializar matrizes do pre-condiionador
Low=[];
Upp=[];

A=gamma_diff*Ad+u_convec*Ac+factor*I;

iter=0;
residuo=1;
numero_condicao=cond(A);
b = b1';
rrrr = A\b1';
    
while residuo>10^-12 && iter<2 % Este ciclo não é estritamente necessário, apenas caso se utilize correcção diferida
    iter=iter+1;
    if iter==1
         rhs=b1'-(-factor*I)*sol_ini;
    else
         rhs=b1'-(-factor*I)*phi_num;
    end
%     aux=inv(gamma_diff*Ad+u_convec*Ac+factor*I)*rhs;
%       aux=inv(gamma_diff*Ad)*rhs;


% Resolve-se directamente o sistema com o solver bicgstab
    [Low,Upp] = ilu(A,struct('type','nofill','droptol',1e-12));
    [aux,fl0,rr0,it0,rv0]=bicgstab(A,rhs,10^-12,10^6,Low,Upp,sol_ini);

    res_max=0;
    res_total=0;
    res_total2=0;
    if iter ~= 1
        for i=1:cell_num
            rescell(i)=abs(aux(i)-phi_num(i));
            if rescell(i)>res_max
                res_max=rescell(i);
            end
            res_total=res_total+rescell(i)*cell_vol(i);
            res_total2=res_total2+(rescell(i)^2*cell_vol(i));
        end
        res(iter,1)=res_total/L;
        res(iter,2)=sqrt(res_total2/L);
        res(iter,3)=res_max;
        residuo=res(iter,2);
    end
    if rem(iter,100)==0
        fprintf('\n%d %E %E %E\n',iter,res(iter,1),res(iter,2),res(iter,3));
    end
    phi_num=aux;
 
end

if residuo>10^-9 || res(iter,3)>0.05
    fl0=5;
end

fprintf('\n numero de iterações=%d res_norma1=%E res_norma2=%E res_norma3=%E\n',iter,res(iter,1),res(iter,2),res(iter,3));
%
end
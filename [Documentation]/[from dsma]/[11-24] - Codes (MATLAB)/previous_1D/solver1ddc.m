clear all
clc
close all

%
%% Declara��o de Variavies %%
%
global u_convec gamma_diff numero_condicao fl0;
global L Lref cell_num face_num;
global x faces cell_vol;
global phi lap_phi phi_West phi_East;
global phi_num lap_phi_num;
global norma1_phi erro_phi_max erro_phi norma1_lap erro_lap_max erro_lap
%
%% Dados de Input %%
%
uniforme=true;                       % Tipo de Malha -----------------------> true----->Malha Uniforme       false---->Malha N�o Uniforme %
metodo1='4thOrder';                  % Metedo Numerico Alta Ordem ----------> 2ndOrder->M�todo de 2� Ordem   4thOrder->M�todo de 4� Ordem    6thOrder->M�todo de 6� Ordem    8thOrder->M�todo de 8� Ordem %
% metodo2='4thOrder';                  % Metedo Numerico Baixa Ordem ---------> 2ndOrder->M�todo de 2� Ordem   4thOrder->M�todo de 4� Ordem    6thOrder->M�todo de 6� Ordem    8thOrder->M�todo de 8� Ordem %
metodo2 = metodo1;

solution='sin';                      % Solu��o Numerica que se est� a usar -> sin------>sin(3*pi*x)          exp------>exp(-(x-mu)^2/s) %
plots=true;                         % Plot dos Resultados %
%
cell_num=10;                        % Numero de C�lulas em que o dominio deve ser dividido %
face_num=cell_num+1;                 % Numero de Faces do Dominio %
%
L=1;                                 % Comprimento do Dominio LxL %
Lref=L/(cell_num);                   % Comprimento de Referencia da Malha %
%
u_convec = 0;                         % Parametro da velocidade convectiva 
gamma_diff = -0.1;                     % Parametro da contribui��o difusiva

Peclet_number = u_convec/gamma_diff*Lref;
%% Prints da Informa��o %%
%
fprintf('Solver 1D com Corre��o Diferida\n');
if uniforme
    fprintf('Malha Uniforme\t');
else
    fprintf('Malha N�o Uniforme\t');
end
%
if metodo1=='2ndOrder'
    fprintf('M�todo de 2� Ordem ');
elseif metodo1=='4thOrder'
    fprintf('M�todo de 4� Ordem ');
elseif metodo1=='6thOrder'
    fprintf('M�todo de 6� Ordem ');
elseif metodo1=='8thOrder'
    fprintf('M�todo de 8� Ordem ');
else
    error('\n\n\tERRO: M�todo Desconhecido\n\n');
end
%
% if metodo2=='2ndOrder'
%     fprintf('em Corre��o Diferida de 2� Ordem\t');
% elseif metodo2=='4thOrder'
%     fprintf('em Corre��o Diferida de 4� Ordem\t');
% elseif metodo2=='6thOrder'
%     fprintf('em Corre��o Diferida de 6� Ordem\t');
% elseif metodo2=='8thOrder'
%     fprintf('em Corre��o Diferida de 8� Ordem\t');
% else
%     error('\n\n\tERRO: M�todo Desconhecido\n\n');
% end
% %
fprintf('Dominio com %d C�lulas\n',cell_num);
%% Gera��o da Malha %%
%
fprintf('\n\nConstru��o da Malha\t\t\t\t\t\t Inicio \t');
[x,faces,cell_vol]=geradormalha(uniforme);
fprintf('... \t Fim\n');
%
%% Solu��o Analitica %%
%
fprintf('\n\nCalculo dos Valores Analiticos\t\t\t Inicio \t');
[phi,lap_phi,phi_West,phi_East]=analyticalsolution(solution,metodo1);
fprintf('... \t Fim\n');
%
%% Corre��o Diferida %%

fprintf('\n\nC�lculo dos Valores Numericos\t\t\t Inicio \t');
[phi_num,iter,res]=DeferredCorrection(metodo1,metodo2);
fprintf('... \t Fim\n');
%
%% Calculo do Erro %%
%
fprintf('\n\nC�lculo do Erro\t\t\t\t\t\t\t Inicio \t');
[norma1_phi,erro_phi_max,erro_phi]=errorcalculationDC;
fprintf('... \t Fim\n');
%
%% Plot dos Resultados
%
if plots
    fprintf('\n\nPlot dos Resultados\t\t\t\t\t\t Inicio \t');
    %
    figure (1)
    plot(x,phi,x,phi_num,'--');
    legend('\phi_{analitico}','\phi_{numerico}');
    title('\phi');
    %
    figure (2)
    plot(x,erro_phi);
    title('\epsilon(\phi)');
    %
    % Plot da Malha %
    y=[0;1];
    [X,Y]=meshgrid(faces,y);
    figure (3)
    plot(X,Y,'b',x,0.5,'kx');
    fprintf('... \t Fim\n');
    %
    figure (4)
    subplot(2,1,1);
    plot(x,phi,x,phi_num,'--');
    legend('\phi_{analitico}','\phi_{numerico}','Location', 'Best');
    title('\phi');
    subplot(2,1,2);
    plot(x,erro_phi);
    legend('\epsilon(\phi)');
    %
    figure (5)
    plot((1:iter),res(:,1),(1:iter),res(:,2),'r--',(1:iter),res(:,3),'k-.');
    title('Residuo');
    legend('||r||_1','||r||_2','||r||_{\infty}');
end
%
%% Escrita dos Resultados para ficheiro %%
% %
% fprintf('\n\nEscrita dos Resultados\t\t\t\t\t Inicio \t');
% filetitle=['results',num2str(cell_num),'.dat'];
% fid=fopen(filetitle,'w');
% fprintf(fid,'\t\t\tn=%d\tLref=%E\tnorma1_phi=%E\tnormamax_phi=%E\n',cell_num,Lref,norma1_phi,erro_phi_max);
% fclose(fid);
% fprintf('... \t Fim\n');
%
%% Apresenta��o de Resultados %%
%
fprintf('\n\nResultados\n');
fprintf('\t\t\tn=%d\tLref=%E\tnorma1_phi=%E\tnormamax_phi=%E\n',cell_num,Lref,norma1_phi,erro_phi_max);
fprintf('\t\t\t%E %E\n',norma1_phi,erro_phi_max);
Peclet_number
fprintf('Fim de Execu��o \n');
%
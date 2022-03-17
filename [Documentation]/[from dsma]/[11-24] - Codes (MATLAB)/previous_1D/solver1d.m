clear all
clc

%
TempoGlobal=cputime;
%
%% Declara��o de Variavies %%
%
global TempoInverterMatriz;
global L Lref cell_num face_num;
global x faces cell_vol;
global phi lap_phi phi_West phi_East flux_phi_West flux_phi_East;
global phi_num lap_phi_num;
global norma1_phi erro_phi_max erro_phi norma1_lap erro_lap_max erro_lap
%
%% Dados de Input %%
%
uniforme=true;                                              % Tipo de Malha -----------------------> true----->Malha Uniforme       false---->Malha N�o Uniforme %
explicito=false;                                            % Tipo de C�lculo ---------------------> true----->C�lculo Explicito    false---->C�lculo Implicito %
metodo='8thOrder';                                          % Metedo Numerico ---------------------> 2ndOrder->M�todo de 2� Ordem   4thOrder->M�todo de 4� Ordem    6thOrder->M�todo de 6� Ordem    8thOrder->M�todo de 8� Ordem %
solution='sin';                                             % Solu��o Numerica que se est� a usar -> sin------>sin(3*pi*x)          exp------>exp(-(x-mu)^2/s) %
dirichlet=false;                                             % Condi��o de Fronteira ---------------> true----->Dirichlet            false---->Neumann % 
plots=true;                                                % Plot dos Resultados %
%
cell_num=400;                                               % Numero de C�lulas em que o dominio deve ser dividido %
face_num=cell_num+1;                                        % Numero de Faces do Dominio %
%
L=1;                                                        % Comprimento do Dominio LxL %
Lref=L/(cell_num);                                          % Comprimento de Referencia da Malha %
%
%% Prints da Informa��o %%
%
fprintf('Solver 1D\n');
if uniforme
    fprintf('Malha Uniforme\t');
else
    fprintf('Malha N�o Uniforme\t');
end
%
if explicito
    fprintf('C�lculo Explicito\t');
else
    fprintf('C�lculo Implicito\t');
end
if dirichlet
    fprintf('Condi��o de Fronteira de Dirichlet\t');
else
    fprintf('Condi��o de Fronteira de Neumann\t');
end
if metodo=='2ndOrder'
    fprintf('M�todo de 2� Ordem\t');
elseif metodo=='4thOrder'
    fprintf('M�todo de 4� Ordem\t');
elseif metodo=='6thOrder'
    fprintf('M�todo de 6� Ordem\t');
elseif metodo=='8thOrder'
    fprintf('M�todo de 8� Ordem\t');
else
    error('\n\n\tERRO: M�todo Desconhecido\n\n');
end
%
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
[phi,lap_phi,phi_West,phi_East,flux_phi_West,flux_phi_East]=analyticalsolution(solution,metodo);
fprintf('... \t Fim\n');
%
%% Solu��o Numerica %%
%
fprintf('\n\nC�lculo dos Valores Numericos\t\t\t Inicio \t');
temposolver=cputime;
if metodo=='2ndOrder'
    [phi_num,lap_phi_num]=Scheme2ndOrder(explicito,dirichlet);
elseif metodo=='4thOrder'
    [phi_num,lap_phi_num]=Scheme4thOrder(explicito,dirichlet);
elseif metodo=='6thOrder'
    [phi_num,lap_phi_num]=Scheme6thOrder(explicito,dirichlet);
elseif metodo=='8thOrder'
    [phi_num,lap_phi_num]=Scheme8thOrder(explicito,dirichlet);
else
    error('\n\n\tERRO: M�todo Desconhecido\n\n');
end
%
temposolver=cputime-temposolver;
fprintf('... \t Fim\n');
%
%% Calculo do Erro %%
%
fprintf('\n\nC�lculo do Erro\t\t\t\t\t\t\t Inicio \t');
[norma1_phi,erro_phi_max,erro_phi,norma1_lap,erro_lap_max,erro_lap,norma1SF,erromaximoSF]=errorcalculation(explicito);
fprintf('... \t Fim\n');

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
    plot(x,lap_phi,x,lap_phi_num,'--');
    title('\nabla^2\phi');
    legend('\nabla^2\phi_{analitico}','\nabla^2\phi_{numerico}');
    %
    figure (3)
    plot(x,erro_phi);
    title('\epsilon(\phi)');
    %
    figure (4)
    plot(x,erro_lap);
    title('\epsilon(\nabla^2\phi)');
    %
    % Plot da Malha %
    y=[0;1];
    [X,Y]=meshgrid(faces,y);
    figure (5)
    plot(X,Y,'b',x,0.5,'x');
    fprintf('... \t Fim\n');
    %
    figure (6)
    subplot(2,1,1);
    plot(x,phi,x,phi_num,'--');
    legend('\phi_{analitico}','\phi_{numerico}');
    title('\phi');
    subplot(2,1,2);
    plot(x,erro_phi);
    legend('\epsilon(\phi)');
    %
    figure (7)
    subplot(2,1,1);
    plot(x,lap_phi,x,lap_phi_num,'--');
    legend('\nabla^2\phi_{analitico}','\nabla^2\phi_{numerico}');
    title('\nabla^2\phi');
    subplot(2,1,2);
    plot(x,erro_lap);
    legend('\epsilon(\nabla^2\phi)');
end
%
%% Escrita dos Resultados para ficheiro %%
%
fprintf('\n\nEscrita dos Resultados\t\t\t\t\t Inicio \t');
filetitle=['results',num2str(cell_num),'.dat'];
fid=fopen(filetitle,'w');
fprintf(fid,'\t\t\tn=%d\tLref=%E\tnorma1_phi=%E\tnormamax_phi=%E\n',cell_num,Lref,norma1_phi,erro_phi_max);
fprintf(fid,'\t\t\tn=%d\tLref=%E\tnorma1_lap=%E\tnormamax_lap=%E\n',cell_num,Lref,norma1_lap,erro_lap_max);
fclose(fid);
fprintf('... \t Fim\n');
%
TempoGlobal=cputime-TempoGlobal;
%
%% Apresenta��o de Resultados %%
%
fprintf('\n\nResultados\n');
fprintf('\t\t\tn=%d\tLref=%E\tnorma1_phi=%E\tnormamax_phi=%E\n',cell_num,Lref,norma1_phi,erro_phi_max);
fprintf('\t\t\tn=%d\tLref=%E\tnorma1_lap=%E\tnormamax_lap=%E\n',cell_num,Lref,norma1_lap,erro_lap_max);
fprintf('\t\t\tTempo Global de Execu��o\t\t%E\n\t\t\tTempo de Resolu��o da Matriz\t%E\n\t\t\tTempo de Execu��o do Solver\t\t%E\n',TempoGlobal,TempoInverterMatriz,temposolver);
fprintf('Fim de Execu��o \n');
%
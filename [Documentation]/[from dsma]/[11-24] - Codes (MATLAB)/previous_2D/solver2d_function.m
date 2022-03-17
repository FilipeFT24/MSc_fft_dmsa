function [Lref, norma1_phi, erro_phi_max, erro_medio_poly, erro_max_poly]=solver2d_function(metodo1)
%
%% Declara��o de Variavies %%
%
global TempoGlobal tempo_malha tempo_anal tempo_num tempo_A tempo_stencil tempo_rec tempo_gmres tempo_erro tempo_plots;
global fid;
global L Lref cell_side vert_side cell_num face_num vert_num;
global verts cells faces;
global cell_verts cell_faces cell_vol cell_norm cell_bound cell_vert_num cell_face_num;
global face_vert face_cells face_area face_bound;
global vert_cells vert_cell_num vert_face_num;
global phi lap_phi  phi_faces flux_phi_faces;
global phi_num lap_phi_num A source source_faces source_cells stencil_cells stencil_faces stencil_size T D;
global norma1_phi erro_phi_max erro_phi norma1_lap erro_lap_max erro_lap X;
global u_convec_x u_convec_y gamma_diff order fig face_w_gauss mix_method G solution T_border restos constrained_source extended_stencil increase_gauss_points exact_coeffs neuman dimensional_correction;
%
TempoGlobal=cputime;
%
%% Dados de Input %%
%
malha='cart';                               % Tipo de Malha ------------> cart-->Cartesiana  tri-->Triangular
solution='sin';                             % Solu��o Numerica que se est� a usar -> sin--->sin(3*pi*x)             exp--->exp(-(x-mu)^2/s) %
equation='diffusion';                       % Equa��o que se pretende resolver
metodo=metodo1;                             % Metedo Numerico ---------------------> FDM_2->2� Ordem                WLS_2->M�todo de 2� Ordem    WLS_4->M�todo de 4� Ordem WLS_6->M�todo de 6� Ordem    WLS_8->M�todo de 8� Ordem %
%
uniforme=true;                              % Tipo de Malha Cartesiana -> treu-->Uniforme    false-->N�o Uniforme
explicito=false;                            % Tipo de C�lculo ---------------------> true-->C�lculo Explicito       false->C�lculo Implicito %
dirichlet=true;                             % Condi��o de Fronteira ---------------> true-->Dirichlet               false->Neumann %
ponderado=true;                             % Pondera��o nos Minimos Quadrados-----> true-->Pondera��o              false->Sem Pondera��o w=1;
GMRES=true;                                 % Solver ------------------------------> true-->GMRES                   false->BICGSTAB %
ILU=true;                                   % Pre Condicionador -------------------> true-->Precondicionador ILU    false->N�o usa Precondicionador %
splots=false;                                 % Plot dos Resultados %
%
% face_w_gauss = true;                        % Nas fronteiras utiliza(true) ou nao(false) os pontos de gauss
% mix_method = true;

L=1;
%
fid=fopen('resultados.txt','w');
%
%cell_side=20;
vert_side=cell_side+1;
%
%% Prints da Informa��o %%
%
informacao(malha,solution,equation,metodo,explicito,uniforme,dirichlet,ponderado,GMRES,ILU);
%
%% Gera��o da Malha %%
%
fprintf('\n\nConstru��o da Malha\t\t\t\t\t\t Inicio ... ');
fprintf(fid,'\n\nConstru��o da Malha\t\t\t\t\t\t Inicio ... ');
tempo_malha=cputime;
%
if strcmp(malha,'cart')==1
    [verts,cells,faces,cell_verts,cell_faces,cell_vol,face_area,face_vert,cell_norm,cell_num,vert_num,face_num,Lref]=CartMesh1(uniforme);
    [face_cells,vert_cells,face_bound,cell_bound,cell_vert_num, cell_face_num,vert_cell_num,vert_face_num]=CartMesh2;
elseif strcmp(malha,'tr')==1
    error('\n\n\tERRO: Fun��o Indisponivel\n\n');
else
    error('\n\n\tERRO: Malha Desconhecida\n\n');
end
%
tempo_malha=cputime-tempo_malha;
tempo_total=cputime-TempoGlobal;
%
fprintf('%f ... %f ',tempo_malha,tempo_total);
fprintf('... Fim\n');
fprintf(fid,'%f ... %f ',tempo_malha,tempo_total);
fprintf(fid,'... Fim\n');
%
%l% Solu��o Analitica %%
%
fprintf('\n\nCalculo dos Valores Analiticos\t\t\t Inicio ... ');
fprintf(fid,'\n\nCalculo dos Valores Analiticos\t\t\t Inicio ... ');
tempo_anal=cputime;
%
[phi,lap_phi,phi_faces,flux_phi_faces]=AnalyticalSolution(solution,metodo,equation);
%
tempo_anal=cputime-tempo_anal;
tempo_total=cputime-TempoGlobal;
%
fprintf('%f ... %f ',tempo_anal,tempo_total);
fprintf('... Fim\n');
fprintf(fid,'%f ... %f ',tempo_anal,tempo_total);
fprintf(fid,'... Fim\n');
%
%% Solu��o Num�rica
%
fprintf('\n\nC�lculo dos Valores Numericos\t\t\t Inicio ... |||||||| ... |||||||| ... |||\n');
fprintf(fid,'\n\nC�lculo dos Valores Numericos\t\t\t Inicio ... |||||||| ... |||||||| ... |||\n');
tempo_num=cputime;
%
if metodo=='FDM_2'
    if strcmp(malha,'cart')==1
        [phi_num,lap_phi_num,A,source,source_faces,tempo_A,tempo_gmres]=FiniteDiferenceMethod2ndOrder(explicito,dirichlet,equation,GMRES,ILU);
    else
        error('\n\nERRO: N�o Implentado para Malhas Triangulares\n\n');
    end
elseif metodo=='WLS_2'
    order=2;
    [phi_num,lap_phi_num,A,source,source_faces,source_cells,stencil_cells,stencil_faces,stencil_size,T,D,tempo_stencil,tempo_rec, tempo_A,tempo_gmres]=WeightedLeastSquares(order,dirichlet,equation,GMRES,ILU,ponderado);
elseif metodo=='CDM_2'
    error('\n\n\tERRO: M�todo N�o Implementado\n\n');
    %[phi_num,lap_phi_num,A,Aw,source,sourcew,source_faces,stencil_cells,stencil_faces,stencil_size,T,D]=DC2ndOrder(dirichlet);
elseif metodo=='WLS_4'
    order=4;
    [phi_num,lap_phi_num,A,source,source_faces,source_cells,stencil_cells,stencil_faces,stencil_size,T,D,tempo_stencil,tempo_rec, tempo_A,tempo_gmres]=WeightedLeastSquares(order,dirichlet,equation,GMRES,ILU,ponderado);
elseif metodo=='WLS_6'
    order=6;
    [phi_num,lap_phi_num,A,source,source_faces,source_cells,stencil_cells,stencil_faces,stencil_size,T,D,tempo_stencil,tempo_rec, tempo_A,tempo_gmres]=WeightedLeastSquares(order,dirichlet,equation,GMRES,ILU,ponderado);
elseif metodo=='WLS_8'
    order=8;
    [phi_num,lap_phi_num,A,source,source_faces,source_cells,stencil_cells,stencil_faces,stencil_size,T,D,tempo_stencil,tempo_rec, tempo_A,tempo_gmres]=WeightedLeastSquares(order,dirichlet,equation,GMRES,ILU,ponderado);
else
    error('\n\n\tERRO: M�todo Desconhecido\n\n');
end
%
tempo_num=cputime-tempo_num;
tempo_total=cputime-TempoGlobal;
%
fprintf('\n\nC�lculo dos Valores Numericos\t\t\t Inicio ... ');
fprintf('%f ... %f ',tempo_num, tempo_total);
fprintf('... Fim\n');
fprintf(fid,'\n\nC�lculo dos Valores Numericos\t\t\t Inicio ... ');
fprintf(fid,'%f ... %f ',tempo_num, tempo_total);
fprintf(fid,'... Fim\n');
%%% Calculo do Erro %%
%
tempo_erro=cputime;
fprintf('\n\nCalculo do Erro\t\t\t\t\t\t\t Inicio ... ');
fprintf(fid,'\n\nCalculo do Erro\t\t\t\t\t\t\t Inicio ... ');
%
[norma1_phi,norma1_lap,erro_phi_max,erro_lap_max,erro_phi,erro_lap,X]=errorcalculation(explicito);
 [erro_medio_poly,erro_max_poly]= Polynomial_total_reconstruction();
% erro_medio_poly=0;
%  erro_max_poly=0;

%
tempo_erro=cputime-tempo_erro;
tempo_total=cputime-TempoGlobal;
%
fprintf('%f ... %f ',tempo_erro, tempo_total);
fprintf('... Fim\n');
fprintf(fid,'%f ... %f ',tempo_erro, tempo_total);
fprintf(fid,'... Fim\n');
%
%% Plot dos Resultados %%
%
if splots
    fprintf('\n\nPlot dos Resultados\t\t\t\t\t\t Inicio ... ');
    fprintf(fid,'\n\nPlot dos Resultados\t\t\t\t\t\t Inicio ... ');
    tempo_plots=cputime;
    %
    plots();
    %
    tempo_plots=cputime-tempo_plots;
    tempo_total=cputime-TempoGlobal;
    %
    fprintf('%f ... %f ',tempo_plots, tempo_total);
    fprintf('... Fim\n');
    fprintf(fid,'%f ... %f ',tempo_plots, tempo_total);
    fprintf(fid,'... Fim\n');
else
    fprintf('\n\nPlot dos Resultados\t\t\t\t\t\t Inicio ...');
    fprintf('N�o ... Fim\n');
    fprintf(fid,'\n\nPlot dos Resultados\t\t\t\t\t\t Inicio ...');
    fprintf(fid,'N�o ... Fim\n');
end
%
%% Export dos Resultados %%
%
fprintf('\n\nExportar Resultados\t\t\t\t\t\t Inicio ... ');
fprintf(fid,'\n\nExportar Resultados\t\t\t\t\t\t Inicio ... ');
tempo_export=cputime;
%
exportplots(explicito);
%
tempo_export=cputime-tempo_export;
TempoGlobal=cputime-TempoGlobal;
%
fprintf('%f ... %f ',tempo_export, TempoGlobal);
fprintf('... Fim\n');
fprintf(fid,'%f ... %f ',tempo_export, TempoGlobal);
fprintf(fid,'... Fim\n');
%
%% Apresenta��o de Resultados %%
%
fprintf('\n\nResultados\n');
fprintf('\t\t\tn=%d\tLref=%E\tnorma1_phi=%E\tnormamax_phi=%E\t%f %f\n',cell_num,Lref,norma1_phi,erro_phi_max,X(2,1),X(2,2));
fprintf('\t\t\tn=%d\tLref=%E\tnorma1_lap=%E\tnormamax_lap=%E\t%f %f\n',cell_num,Lref,norma1_lap,erro_lap_max,X(1,1),X(1,2));
fprintf('\t\t\tTempo Global de Execu��o\t\t%E\n',TempoGlobal);
fprintf('Fim de Execu��o \n');
%
fprintf(fid,'\n\nResultados\n');
fprintf(fid,'\t\t\tn=%d\tLref=%E\tnorma1_phi=%E\tnormamax_phi=%E\t%f %f\n',cell_num,Lref,norma1_phi,erro_phi_max,X(2,1),X(2,2));
fprintf(fid,'\t\t\tn=%d\tLref=%E\tnorma1_lap=%E\tnormamax_lap=%E\t%f %f\n',cell_num,Lref,norma1_lap,erro_lap_max,X(1,1),X(1,2));
fprintf(fid,'\t\t\tTempo Global de Execu��o\t\t%E\n',TempoGlobal);
fprintf(fid,'Fim de Execu��o \n');
%
fclose(fid);
end
function informacao(malha,solution,equation,metodo,explicito,uniforme,dirichlet,ponderado,GMRES,ILU)
%
%% Declara��o de Variavies %%
%
global L Lref cell_side vert_side cell_num face_num vert_num;
global verts cells faces;
global cell_verts cell_faces cell_vol cell_norm cell_bound cell_vert_num cell_face_num;
global face_vert face_cells face_area face_bound;
global vert_cells vert_cell_num vert_face_num;
global phi lap_phi  phi_faces flux_phi_faces;
global phi_num lap_phi_num A source source_faces stencil_cells stencil_faces stencil_size T D;
global norma1_phi erro_phi_max erro_phi norma1_lap erro_lap_max erro_lap X;
global TempoGlobal tempo_malha tempo_anal tempo_num tempo_A tempo_stencil tempo_rec tempo_gmres tempo_erro tempo_plots;
global fid;
%
%% Prints da Informa��o %%
%
fprintf('Solver 2D\n');
fprintf(fid,'Solver 2D\n');
%
if strcmp(equation,'diffusion')==1
    fprintf('Equa��o Difusiva\t\t');
    fprintf(fid,'Equa��o Difusiva\t\t');
elseif strcmp(equation,'convection')==1
    fprintf('Equa��o Convectiva\t\t');
    fprintf(fid,'Equa��o Convectiva\t\t');
elseif strcmp(equation,'reaction')==1
    fprintf('Equa��o Reac��o Difus�o\t\t');
    fprintf(fid,'Equa��o Reac��o Difus�o\t\t');
elseif strcmp(equation,'convdif')==1
    fprintf('Equa��o Convec��o Difus�o\t');
    fprintf(fid,'Equa��o Convec��o Difus�o\t');
else 
    error('\n\n\tERRO: Equa��o Desconhecida\n\n');
end
%
if metodo=='FDM_2'
    fprintf('M�todo de Diferen�as Finitas 2� Ordem ');
    fprintf(fid,'M�todo de Diferen�as Finitas 2� Ordem ');
elseif metodo=='WLS_2'
    fprintf('M�todo Minimos Quadrados 2� Ordem ');
    fprintf(fid,'M�todo Minimos Quadrados 2� Ordem ');
elseif metodo=='CDM_2'
    fprintf('M�todo Corre��o Diferida 2� Ordem ');
    fprintf(fid,'M�todo Corre��o Diferida 2� Ordem ');
elseif metodo=='WLS_4'
    fprintf('M�todo Minimos Quadrados 4� Ordem ');
    fprintf(fid,'M�todo Minimos Quadrados 4� Ordem ');
elseif metodo=='WLS_6'
    fprintf('M�todo Minimos Quadrados 6� Ordem ');
    fprintf(fid,'M�todo Minimos Quadrados 6� Ordem ');
elseif metodo=='WLS_8'
    fprintf('M�todo Minimos Quadrados 8� Ordem ');
    fprintf(fid,'M�todo Minimos Quadrados 8� Ordem ');
else
    error('\n\n\tERRO: M�todo Desconhecido\n\n');
end
%
if strcmp(metodo,'FDM_2')==0
    if ponderado==true
        fprintf('Ponderado ');
        fprintf(fid,'Ponderado ');
    else
        fprintf('N�o Ponderado ');
        fprintf(fid,'N�o Ponderado ');
    end
end
%
if explicito
    fprintf('Explicito\n');
    fprintf(fid,'Explicito\n');
else
    fprintf('Implicito\n');
    fprintf(fid,'Implicito\n');
end
%
if strcmp(malha,'cart')==1
    fprintf('Malha Cartesiana\t\t');
    fprintf(fid,'Malha Cartesiana\t\t');
elseif strcmp(malha,'tri')==1
    fprintf('Malha Triangular\t\t');
    fprintf(fid,'Malha Triangular\t\t');
else
    error('\n\n\tERRO: Malha Desconhecida\n\n');
end
%
if dirichlet
    fprintf('Condi��es de Fronteira de Dirichlet\n');
    fprintf(fid,'Condi��es de Fronteira de Dirichlet\n');
else
    fprintf('Condi��es de Fronteira de Neumann\n');
    fprintf(fid,'Condi��es de Fronteira de Neumann\n');
end
%
if strcmp(malha,'cart')==1
    if uniforme
        fprintf('Malha Uniforme\t\t\t');
        fprintf(fid,'Malha Uniforme\t\t\t');
    else
        fprintf('Malha N�o Uniforme\t\t');
        fprintf(fid,'Malha N�o Uniforme\t\t');
    end
end
%
fprintf('Dominio com %dx%d C�lulas\n',cell_side,cell_side);
fprintf(fid,'Dominio com %dx%d C�lulas\n',cell_side,cell_side);
%
if strcmp(solution,'sin')==1
    fprintf('Solu��o Sinusoiodal\t\t');
    fprintf(fid,'Solu��o Sinusoiodal\t\t');
elseif strcmp(solution,'exp')==1
    fprintf('Solu��o Exponencial\t\t');
    fprintf(fid,'Equa��o Exponencial\t\t');
elseif strcmp(solution,'2nd')==1
    fprintf('Polinomio Linear\t\t');
    fprintf(fid,'Polinomio Linear\t\t');
elseif strcmp(solution,'4th')==1
    fprintf('Polinomio Cubico\t\t');
    fprintf(fid,'Polinomio Cubico\t\t');
elseif strcmp(solution,'6th')==1
    fprintf('Polinomio 5� Grau\t\t');
    fprintf(fid,'Polinomio 5� Grau\t\t');
elseif strcmp(solution,'8th')==1
    fprintf('Polinomio 7� Grau\t\t');
    fprintf(fid,'Polinomio 7� Grau\t\t');
else
    fprintf('Outra Solu��o\t\t');
    fprintf(fid,'Outra Solu��o\t\t');
end
%
if GMRES
    fprintf('Solver GMRES\t\t');
    fprintf(fid,'Solver GMRES\t\t');
else
    fprintf('Solver BICGSTAB\t\t');
    fprintf(fid,'Solver BICGSTAB\t\t');
end
%
if ILU
    fprintf('Precondicionador ILU\n');
    fprintf(fid,'Precondicionador ILU\n');
else
    fprintf('Sem Precondicionador\n');
    fprintf(fid,'Sem Precondicionador\n');
end
%
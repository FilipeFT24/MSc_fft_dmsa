function informacao(malha,solution,equation,metodo,explicito,uniforme,dirichlet,ponderado,GMRES,ILU)
%
%% Declaração de Variavies %%
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
%% Prints da Informação %%
%
fprintf('Solver 2D\n');
fprintf(fid,'Solver 2D\n');
%
if strcmp(equation,'diffusion')==1
    fprintf('Equação Difusiva\t\t');
    fprintf(fid,'Equação Difusiva\t\t');
elseif strcmp(equation,'convection')==1
    fprintf('Equação Convectiva\t\t');
    fprintf(fid,'Equação Convectiva\t\t');
elseif strcmp(equation,'reaction')==1
    fprintf('Equação Reacção Difusão\t\t');
    fprintf(fid,'Equação Reacção Difusão\t\t');
elseif strcmp(equation,'convdif')==1
    fprintf('Equação Convecção Difusão\t');
    fprintf(fid,'Equação Convecção Difusão\t');
else 
    error('\n\n\tERRO: Equação Desconhecida\n\n');
end
%
if metodo=='FDM_2'
    fprintf('Método de Diferenças Finitas 2ª Ordem ');
    fprintf(fid,'Método de Diferenças Finitas 2ª Ordem ');
elseif metodo=='WLS_2'
    fprintf('Método Minimos Quadrados 2ª Ordem ');
    fprintf(fid,'Método Minimos Quadrados 2ª Ordem ');
elseif metodo=='CDM_2'
    fprintf('Método Correção Diferida 2ª Ordem ');
    fprintf(fid,'Método Correção Diferida 2ª Ordem ');
elseif metodo=='WLS_4'
    fprintf('Método Minimos Quadrados 4ª Ordem ');
    fprintf(fid,'Método Minimos Quadrados 4ª Ordem ');
elseif metodo=='WLS_6'
    fprintf('Método Minimos Quadrados 6ª Ordem ');
    fprintf(fid,'Método Minimos Quadrados 6ª Ordem ');
elseif metodo=='WLS_8'
    fprintf('Método Minimos Quadrados 8ª Ordem ');
    fprintf(fid,'Método Minimos Quadrados 8ª Ordem ');
else
    error('\n\n\tERRO: Método Desconhecido\n\n');
end
%
if strcmp(metodo,'FDM_2')==0
    if ponderado==true
        fprintf('Ponderado ');
        fprintf(fid,'Ponderado ');
    else
        fprintf('Não Ponderado ');
        fprintf(fid,'Não Ponderado ');
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
    fprintf('Condições de Fronteira de Dirichlet\n');
    fprintf(fid,'Condições de Fronteira de Dirichlet\n');
else
    fprintf('Condições de Fronteira de Neumann\n');
    fprintf(fid,'Condições de Fronteira de Neumann\n');
end
%
if strcmp(malha,'cart')==1
    if uniforme
        fprintf('Malha Uniforme\t\t\t');
        fprintf(fid,'Malha Uniforme\t\t\t');
    else
        fprintf('Malha Não Uniforme\t\t');
        fprintf(fid,'Malha Não Uniforme\t\t');
    end
end
%
fprintf('Dominio com %dx%d Células\n',cell_side,cell_side);
fprintf(fid,'Dominio com %dx%d Células\n',cell_side,cell_side);
%
if strcmp(solution,'sin')==1
    fprintf('Solução Sinusoiodal\t\t');
    fprintf(fid,'Solução Sinusoiodal\t\t');
elseif strcmp(solution,'exp')==1
    fprintf('Solução Exponencial\t\t');
    fprintf(fid,'Equação Exponencial\t\t');
elseif strcmp(solution,'2nd')==1
    fprintf('Polinomio Linear\t\t');
    fprintf(fid,'Polinomio Linear\t\t');
elseif strcmp(solution,'4th')==1
    fprintf('Polinomio Cubico\t\t');
    fprintf(fid,'Polinomio Cubico\t\t');
elseif strcmp(solution,'6th')==1
    fprintf('Polinomio 5º Grau\t\t');
    fprintf(fid,'Polinomio 5º Grau\t\t');
elseif strcmp(solution,'8th')==1
    fprintf('Polinomio 7º Grau\t\t');
    fprintf(fid,'Polinomio 7º Grau\t\t');
else
    fprintf('Outra Solução\t\t');
    fprintf(fid,'Outra Solução\t\t');
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
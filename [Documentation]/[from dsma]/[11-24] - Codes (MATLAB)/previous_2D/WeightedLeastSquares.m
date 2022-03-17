function [phi_num,lap_phi_num,A,source,source_faces,source_cells,stencil_cells,stencil_faces,stencil_size,T,D,tempo_stencil,tempo_rec, tempo_A,tempo_gmres]=WeightedLeastSquares(order,dirichlet,equation,GMRES,ILU,ponderado)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   26 de Setembro de 2016                                          %
%                                                                                                   %
% Função que implementa a Solução Numérica a partir do Minimos Quadrados de 2ª Ordem                %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global TempoGlobal fid;
global L Lref cell_side vert_side cell_num face_num vert_num;
global verts cells faces;
global cell_verts cell_faces cell_vol cell_norm cell_bound cell_vert_num cell_face_num;
global face_vert face_cells face_area face_bound;
global vert_cells vert_cell_num vert_face_num;
global phi lap_phi  phi_faces flux_phi_faces stencil_cells stencil_faces stencil_size order;
global u_convec_x u_convec_y gamma_diff G T_border constrained_source mix_method neuman face_w_gauss increase_gauss_points exact_coeffs dimensional_correction robin;
%
%% Construção do Stencil %%
%
fprintf('\n\n\tConstrução do Stencil\t\t\t\t Inicio ... ');
fprintf(fid,'\n\n\tConstrução do Stencil\t\t\t\t Inicio ... ');
tempo_stencil=cputime;
%
[stencil_cells,stencil_faces,stencil_size]=Stencil(order);
%
tempo_stencil=cputime-tempo_stencil;
tempo_total=cputime-TempoGlobal;
%
fprintf('%f ... %f ',tempo_stencil, tempo_total);
fprintf('... Fim\n');
fprintf(fid,'%f ... %f ',tempo_stencil, tempo_total);
fprintf(fid,'... Fim\n');
%
%% Pontos de Gauss %%
%
fprintf('\n\n\tPontos de Gauss na Face\t\t\t\t Inicio ... ');
fprintf(fid,'\n\n\tPontos de Gauss na Face\t\t\t\t Inicio ... ');
tempo_gauss=cputime;
%

[G]=GaussFace(order); % [G]=GaussFace(ordem) %
%
tempo_gauss=cputime-tempo_gauss;
tempo_total=cputime-TempoGlobal;
%
fprintf('%f ... %f ',tempo_gauss, tempo_total);
fprintf('... Fim\n');
fprintf(fid,'%f ... %f ',tempo_gauss, tempo_total);
fprintf(fid,'... Fim\n');
%
%% Reconstrução do Polinomio  %%
%
fprintf('\n\n\tReconstrução do Polinomio na Face\t Inicio ... ');
fprintf(fid,'\n\n\tReconstrução do Polinomio na Face\t Inicio ... ');
tempo_rec=cputime;
%
if order==2
    [T,D]=Reconstruction2ndOrder(stencil_cells,stencil_faces,stencil_size,ponderado);
    if face_w_gauss
     [T_border,D_border, constrained_source]=Reconstruction2ndOrder_borders_v2(stencil_cells,stencil_faces,stencil_size,ponderado);
    end
elseif order==4
    [T,D]=Reconstruction4thOrder(stencil_cells,stencil_faces,stencil_size,ponderado);
    if face_w_gauss
    [T_border,D_border, constrained_source]=Reconstruction4thOrder_borders_v2(stencil_cells,stencil_faces,stencil_size,ponderado);
    end
elseif order==6
    if exact_coeffs
    stencil_adpat_excat();
    [T,D]=Reconstruction6thOrder_direct(stencil_cells,stencil_faces,stencil_size,ponderado);
    elseif face_w_gauss
    [T,D]=Reconstruction6thOrder(stencil_cells,stencil_faces,stencil_size,ponderado);
    [T_border,D_border, constrained_source]=Reconstruction6thOrder_borders_v2(stencil_cells,stencil_faces,stencil_size,ponderado);
    else
    [T,D]=Reconstruction6thOrder(stencil_cells,stencil_faces,stencil_size,ponderado);
    end
elseif order==8
    if exact_coeffs
    stencil_adpat_excat();
    [T,D]=Reconstruction8thOrder_direct(stencil_cells,stencil_faces,stencil_size,ponderado);
    elseif face_w_gauss
    [T,D]=Reconstruction8thOrder(stencil_cells,stencil_faces,stencil_size,ponderado);
    [T_border,D_border, constrained_source]=Reconstruction8thOrder_borders_v2(stencil_cells,stencil_faces,stencil_size,ponderado);
    else
        [T,D]=Reconstruction8thOrder(stencil_cells,stencil_faces,stencil_size,ponderado);
    end
else 
    error('\n\nERRO: Reconstrução Não Implementada\n\n');
end
%

pause(1)
tempo_rec=cputime-tempo_rec;
tempo_total=cputime-TempoGlobal;
%
fprintf('%f ... %f ',tempo_rec, tempo_total);
fprintf('... Fim\n');
fprintf(fid,'%f ... %f ',tempo_rec, tempo_total);
fprintf(fid,'... Fim\n');
%
%% Pontos de Gauss %%
if increase_gauss_points && order==4
[G]=GaussFace(order+2); % [G]=GaussFace(ordem) %
end
%% Construção da Matriz A %%
%
fprintf('\n\n\tContrução da Matriz A\t\t\t\t Inicio ...');
tempo_A=cputime;
%
if strcmp(equation,'diffusion')==1
    [Adif,sourcedif,source_faces,source_cells]=MatrixDiffusion(dirichlet,order,T,T_border, constrained_source,G,stencil_cells,stencil_faces,stencil_size);
    [Aconv,sourceconv,source_facesconv,source_cells]=MatrixConvection(dirichlet,order,T,T_border,  constrained_source,G,stencil_cells,stencil_faces,stencil_size);
elseif strcmp(equation,'convection')==1
    error('\n\nERRO: Equação Convectiva Não Implementada\n\n');
elseif strcmp(equation,'reaction')==1
    error('\n\nERRO: Equação Difusão Reacção Não Implementada\n\n');
elseif strcmp(equation,'convdif')==1
    error('\n\nERRO: Equação Convecção Difusão Não Implementada\n\n');
end
%
tempo_A=cputime-tempo_A;
tempo_total=cputime-TempoGlobal;
%
fprintf('\n\n\tConstrução da Matriz A\t\t\t\t Inicio ... ');
fprintf('%f ... %f ',tempo_A, tempo_total);
fprintf('... Fim\n');
fprintf(fid,'\n\n\tConstrução da Matriz A\t\t\t\t Inicio ... ');
fprintf(fid,'%f ... %f ',tempo_A, tempo_total);
fprintf(fid,'... Fim\n');
%

    A = gamma_diff*Adif+Aconv;
    source = gamma_diff*sourcedif+sourceconv+source_cells;
%% Solver %%
if GMRES
    %
    % GMRES %
    %
    fprintf('\n\n\tSolver GMRES\t\t\t\t\t\t ');
    tempo_gmres=cputime;
    %
    % Cálculo do PHI atraves do solver GMRES %
    % ver help gmres %
    % http://www.mathworks.com/help/matlab/ref/gmres.html %
    %
    

    if ILU
        % Precondicionador ILU %
        LL=[];
        UU=[];
         [LL,UU] = ilu(A,struct('type','ilutp','droptol',1e-12));
        %
%        phi_num=gmres(A,source',5,1e-12,10000,LL,UU,phi');
        phi_num=gmres(A,source',5,1e-15,10000,LL,UU);
%         phi_num=bicg(A,source',1e-12,10000,LL,UU,phi');
        %
    else
        % Sem Precondionador %
        %
        phi_num=gmres(A,source',5,1e-10,size(A,1));
        %
    end
    %
    tempo_gmres=cputime-tempo_gmres;
    tempo_total=cputime-TempoGlobal;
    %
    fprintf('\n\n\tSolver GMRES\t\t\t\t\t\t Inicio ... ');
    fprintf('%f ... %f ',tempo_gmres, tempo_total);
    fprintf('... Fim\n');
    fprintf(fid,'\n\n\tSolver GMRES\t\t\t\t\t\t Inicio ... ');
    fprintf(fid,'%f ... %f ',tempo_gmres, tempo_total);
    fprintf(fid,'... Fim\n');
else
    %
    % BICGSTAB %
    %
    fprintf('\n\n\tSolver BICGSTAB\t\t\t\t\t\t ');
    tempo_gmres=cputime;
    %
    % Cálculo do PHI atraves do solver BICGSTAB %
    % ver help bicg %
    % http://www.mathworks.com/help/matlab/ref/bicg.html %
    %
    if ILU
        % Precondicionador ILU %
        %
        [LL,UU] = ilu(A,struct('type','ilutp','droptol',1e-12));
        %
        phi_num=bicg(A,source',1e-15,10000,LL,UU);
    else
        % Sem Precondiocinador %
        phi_num=bicg(A,source',1e-10,size(A,1));
        %
    end
    %
    tempo_gmres=cputime-tempo_gmres;
    tempo_total=cputime-TempoGlobal;
    %
    fprintf('\n\n\tSolver BICSGTAB\t\t\t\t\t\t Inicio ... ');
    fprintf('%f ... %f ',tempo_gmres, tempo_total);
    fprintf('... Fim\n');
    fprintf(fid,'\n\n\tSolver BICSGTAB\t\t\t\t\t\t Inicio ... ');
    fprintf(fid,'%f ... %f ',tempo_gmres, tempo_total);
    fprintf(fid,'... Fim\n');
end
%
% Cálculo do Laplaciano de PHI % TER CUIDADO!!!
lap_phi_num1=A*phi_num+source_faces';
%
for i=1:cell_num
    lap_phi_num(i)=lap_phi_num1(i)/cell_vol(i);
end
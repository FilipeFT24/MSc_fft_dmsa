function [phi_num,lap_phi_num,A,source,source_faces,stencil_cells,stencil_faces,stencil_size,T,D,tempo_stencil,tempo_rec, tempo_A,tempo_gmres]=WLS2ndOrder(dirichlet,equation,GMRES,ILU,ponderado)
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
global phi lap_phi  phi_faces flux_phi_faces;
%
order=2;
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
%% Reconstrução do Polinomio  %%
%
fprintf('\n\n\tReconstrução do Polinomio na Face\t Inicio ... ');
fprintf(fid,'\n\n\tReconstrução do Polinomio na Face\t Inicio ... ');
tempo_rec=cputime;
%
[T,D]=Reconstruction2ndOrder(stencil_cells,stencil_faces,stencil_size,ponderado);
%
tempo_rec=cputime-tempo_rec;
tempo_total=cputime-TempoGlobal;
%
fprintf('%f ... %f ',tempo_rec, tempo_total);
fprintf('... Fim\n');
fprintf(fid,'%f ... %f ',tempo_rec, tempo_total);
fprintf(fid,'... Fim\n');
%
%% Construção da Matriz A %%
%
fprintf('\n\n\tContrução da Matriz A\t\t\t\t Inicio ...');
tempo_A=cputime;
%
[A,source,source_faces]=MatrixDiffWLS(dirichlet,order,T,stencil_cells,stencil_faces,stencil_size);
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
        %
        [LL,UU] = ilu(A,struct('type','ilutp','droptol',1e-5));
        %
        phi_num=gmres(A,source',5,1e-10,size(A,1),LL,UU);
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
        [LL,UU] = ilu(A,struct('type','ilutp','droptol',1e-5));
        %
        phi_num=bicg(A,source',1e-10,size(A,1),LL,UU);
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
% Cálculo do Laplaciano de PHI %
lap_phi_num1=A*phi_num+source_faces';
%
for i=1:cell_num
    lap_phi_num(i)=lap_phi_num1(i)/cell_vol(i);
end
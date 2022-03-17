function [phi_num,lap_phi_num,A,source,source_faces,tempo_A,tempo_gmres]=FiniteDiferenceMethod2ndOrder(explicito,dirichlet,equation,GMRES,ILU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   04 de Julho de 2016                                             %
%                                                                                                   %
% Função que implementa a Solução Numérica a partir do Método de Diferenças Finitas de 2ª Ordem     %
% O Método Explicito ainda não foi implementado, sendo que o método implicito recorre à inversão da % 
% matriz                                                                                            % 
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global TempoGlobal;
global fid;
global L Lref cell_side vert_side cell_num face_num vert_num;
global verts cells faces;
global cell_verts cell_faces cell_vol cell_norm cell_bound cell_vert_num cell_face_num;
global face_vert face_cells face_area face_bound;
global vert_cells vert_cell_num vert_face_num;
global phi lap_phi  phi_faces flux_phi_faces;
%
if explicito
    error('\n\n\tERRO: Método Não Implementado\n\n');
else
    %
    tempo_A=cputime;
    %
    if strcmp(equation,'diffusion')==1
        [A,source,source_faces]=MatrixDiffFDM(dirichlet);
    else
        error('\n\n\tERRO: Método Não Implementado\n\n');
    end
    %
    tempo_A=cputime-tempo_A;
    tempo_total=cputime-TempoGlobal;
    %
    fprintf('\n\n\t\tConstrução da Matriz A\t\t\t Inicio ... ');
    fprintf('%f ... %f ',tempo_A, tempo_total);
    fprintf('... Fim\n');
    fprintf(fid,'\n\n\t\tConstrução da Matriz A\t\t\t Inicio ... ');
    fprintf(fid,'%f ... %f ',tempo_A, tempo_total);
    fprintf(fid,'... Fim\n');
    %
    %% Solver %%
    %
    if GMRES
        %
        % GMRES %
        %
        fprintf('\n\n\t\tSolver GMRES\t\t\t\t\t ');
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
        fprintf('\n\n\t\tSolver GMRES\t\t\t\t\t Inicio ... ');
        fprintf('%f ... %f ',tempo_gmres, tempo_total);
        fprintf('... Fim\n');
        fprintf(fid,'\n\n\t\tSolver GMRES\t\t\t\t\t Inicio ... ');
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
end
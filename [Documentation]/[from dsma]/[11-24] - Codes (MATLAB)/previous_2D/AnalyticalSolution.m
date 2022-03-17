function [phi,lap_phi,phi_faces,flux_phi_faces]=AnalyticalSolution(solution,metodo,equation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Artur Guilherme Vasconcelos                                         %
%                                   28 de Junho de 2016                                             %
%                                                                                                   %
% Função que determina as soluções analiticas                                                       %
%                                                                                                   %
% phi               - Valor da Propriedade no centroide da célula;                                  %
% lap_phi           - Laplaciano da Propriedade na célula, utiliza pontos de Gauss;                 %
% phi_faces         - Valor da Propriedade nas fronteiras;                                          %
% flux_phi_faces    - Valor do Fluxo nas fronteiras, já está a multiplicar pela area da face;       %
%                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global L Lref cell_side vert_side cell_num face_num vert_num;
global verts cells faces;
global cell_verts cell_faces cell_vol cell_norm cell_bound cell_vert_num cell_face_num;
global face_vert face_cells face_area face_bound;
global vert_cells vert_cell_num vert_face_num;
global u_convec_x u_convec_y gamma_diff;
%
%% Celulas do Dominio %%
%
for i=1:cell_num
    %
    % Função Analitica %
    if strcmp(equation,'diffusion')==1
        phi(i)=SolutionDiffusion(solution,cells(i,1),cells(i,2),'anal');
    else
        error('\n\nERRO:Não Implementado\n\n');
    end
    %
    % Calcula o Laplaciano de Phi com base nos Pontos de Gauss %
    if strcmp(equation,'diffusion')==1
        %
        % Equação Difusiva %
        if metodo=='FDM_2' | metodo=='WLS_2'
            %
            % 2ª Ordem %
            if strcmp(solution,'paper')==1
            lap_phi(i)=SolutionDiffusion(solution,cells(i,1),cells(i,2),'lap');   
            else
            lap_phi(i)=gamma_diff*SolutionDiffusion(solution,cells(i,1),cells(i,2),'lap')+u_convec_x*SolutionDiffusion(solution,cells(i,1),cells(i,2),'xflux')+u_convec_y*SolutionDiffusion(solution,cells(i,1),cells(i,2),'yflux');
            end
            %
        elseif metodo=='WLS_4'
            %
            % 4ª Ordem %
            % Calula o Laplaciano de Phi com base nos Pontos de Gauss %
            lap_phi(i)=0;
            %
            % Cada Celula é dividida em Triangulos %
            for j=1:cell_vert_num(i)
                aux1=0;
                if j==cell_vert_num(i)
                    x1=verts(cell_verts(i,j),1);
                    y1=verts(cell_verts(i,j),2);
                    x2=verts(cell_verts(i,1),1);
                    y2=verts(cell_verts(i,1),2);
                    x3=cells(i,1);
                    y3=cells(i,2);
                else
                    x1=verts(cell_verts(i,j),1);
                    y1=verts(cell_verts(i,j),2);
                    x2=verts(cell_verts(i,j+1),1);
                    y2=verts(cell_verts(i,j+1),2);
                    x3=cells(i,1);
                    y3=cells(i,2);
                end
                %
                vol_tri=abs(x1*y2+x3*y1+x2*y3-(x3*y2+x1*y3+x2*y1))/2;
                %
                % Pontos de Gauss %
                [x_gauss,ngauss]=gausspoints(x1,y1,x2,y2,x3,y3,'2D',4);
                %
                for k=1:ngauss
                    if strcmp(solution,'paper')==1
                    aux2=SolutionDiffusion(solution,x_gauss(k,1),x_gauss(k,2),'lap');
                    else
                    aux2=gamma_diff*SolutionDiffusion(solution,x_gauss(k,1),x_gauss(k,2),'lap')+u_convec_x*SolutionDiffusion(solution,x_gauss(k,1),x_gauss(k,2),'xflux')+u_convec_y*SolutionDiffusion(solution,x_gauss(k,1),x_gauss(k,2),'yflux'); 
                    end
                    aux1=aux1+aux2*x_gauss(k,3)*vol_tri/cell_vol(i);
                end
                lap_phi(i)=lap_phi(i)+aux1;
            end
        elseif metodo=='WLS_6'
            %
            % 6ª Ordem %
            %
            % Calula o Laplaciano de Phi com base nos Pontos de Gauss %
            lap_phi(i)=0;
            for j=1:cell_vert_num(i)
                aux1=0;
                if j==cell_vert_num(i)
                    x1=verts(cell_verts(i,j),1);
                    y1=verts(cell_verts(i,j),2);
                    x2=verts(cell_verts(i,1),1);
                    y2=verts(cell_verts(i,1),2);
                    x3=cells(i,1);
                    y3=cells(i,2);
                else
                    x1=verts(cell_verts(i,j),1);
                    y1=verts(cell_verts(i,j),2);
                    x2=verts(cell_verts(i,j+1),1);
                    y2=verts(cell_verts(i,j+1),2);
                    x3=cells(i,1);
                    y3=cells(i,2);
                end
                %
                vol_tri=abs(x1*y2+x3*y1+x2*y3-(x3*y2+x1*y3+x2*y1))/2;
                %
                % Pontos de Gauss %
                [x_gauss,ngauss]=gausspoints(x1,y1,x2,y2,x3,y3,'2D',6);
                %
                for k=1:ngauss
                    if strcmp(solution,'paper')==1
                    aux2=SolutionDiffusion(solution,x_gauss(k,1),x_gauss(k,2),'lap');
                    else
                    aux2=gamma_diff*SolutionDiffusion(solution,x_gauss(k,1),x_gauss(k,2),'lap')+u_convec_x*SolutionDiffusion(solution,x_gauss(k,1),x_gauss(k,2),'xflux')+u_convec_y*SolutionDiffusion(solution,x_gauss(k,1),x_gauss(k,2),'yflux'); 
                    end
                    aux1=aux1+aux2*x_gauss(k,3)*vol_tri/cell_vol(i);
                end
                lap_phi(i)=lap_phi(i)+aux1;
            end
        elseif metodo=='WLS_8'
            %
            % 8ª Ordem %
            %
            % Calula o Laplaciano de Phi com base nos Pontos de Gauss %
            lap_phi(i)=0;
            for j=1:cell_vert_num(i)
                aux1=0;
                if j==cell_vert_num(i)
                    x1=verts(cell_verts(i,j),1);
                    y1=verts(cell_verts(i,j),2);
                    x2=verts(cell_verts(i,1),1);
                    y2=verts(cell_verts(i,1),2);
                    x3=cells(i,1);
                    y3=cells(i,2);
                else
                    x1=verts(cell_verts(i,j),1);
                    y1=verts(cell_verts(i,j),2);
                    x2=verts(cell_verts(i,j+1),1);
                    y2=verts(cell_verts(i,j+1),2);
                    x3=cells(i,1);
                    y3=cells(i,2);
                end
                %
                vol_tri=abs(x1*y2+x3*y1+x2*y3-(x3*y2+x1*y3+x2*y1))/2;
                %
                % Pontos de Gauss %
                [x_gauss,ngauss]=gausspoints(x1,y1,x2,y2,x3,y3,'2D',8);
                %
                for k=1:ngauss                     
                    if strcmp(solution,'paper')==1
                    aux2=SolutionDiffusion(solution,x_gauss(k,1),x_gauss(k,2),'lap');
                    else
                    aux2=gamma_diff*SolutionDiffusion(solution,x_gauss(k,1),x_gauss(k,2),'lap')+u_convec_x*SolutionDiffusion(solution,x_gauss(k,1),x_gauss(k,2),'xflux')+u_convec_y*SolutionDiffusion(solution,x_gauss(k,1),x_gauss(k,2),'yflux'); 
                    end
                    aux1=aux1+aux2*x_gauss(k,3)*vol_tri/cell_vol(i);
                end
                lap_phi(i)=lap_phi(i)+aux1;
            end
        else
            error('\n\nERRO: Não Implementado\n\n');
        end
    else
        error('\n\nERRO: Não Implementado');
    end
end
%
%% Faces de Fronteira %%
%
for i=1:face_num
    if face_bound(i,1)==1
        %
        % Função Analitica %
        if strcmp(equation,'diffusion')==1
            %
            % Equação Difusiva %
            phi_faces(i)=SolutionDiffusion(solution,faces(i,1),faces(i,2),'anal');
        else
            error('\n\nERRO:Não Implementado\n\n');
        end
        %
        % Calula o Fluxo de Phi com base nos Pontos de Gauss %
        if strcmp(equation,'diffusion')==1
            %
            % Equação Difusiva %
            if metodo=='FDM_2' | metodo=='WLS_2'
                %
                % 2ª Ordem %
                aux_x=SolutionDiffusion(solution,faces(i,1),faces(i,2),'xflux');
                aux_y=SolutionDiffusion(solution,faces(i,1),faces(i,2),'yflux');
                %
                flux_phi_faces(i,1)=aux_x;%*face_area(i,2);
                flux_phi_faces(i,2)=aux_y;%*face_area(i,3);
                %
            elseif metodo=='WLS_4'
                %
                % 4ª Ordem %
                x1=verts(face_vert(i,1),1);
                y1=verts(face_vert(i,1),2);
                x2=verts(face_vert(i,2),1);
                y2=verts(face_vert(i,2),2);
                x3=faces(i,1);
                y3=faces(i,2);
                %
                % Pontos de Gauss %
                [x_gauss,ngauss]=gausspoints(x1,y1,x2,y2,x3,y3,'1D',4);
                %
                flux_phi_faces(i,1)=0;
                flux_phi_faces(i,2)=0;
                %
                for j=1:ngauss
                    aux_x=SolutionDiffusion(solution,x_gauss(j,1),x_gauss(j,2),'xflux')*x_gauss(j,3);
                    aux_y=SolutionDiffusion(solution,x_gauss(j,1),x_gauss(j,2),'yflux')*x_gauss(j,3);
                    %
                    flux_phi_faces(i,1)=flux_phi_faces(i,1)+aux_x;
                    flux_phi_faces(i,2)=flux_phi_faces(i,2)+aux_y;
                end
            elseif metodo=='WLS_6'
                %
                % 6ª Ordem %
                x1=verts(face_vert(i,1),1);
                y1=verts(face_vert(i,1),2);
                x2=verts(face_vert(i,2),1);
                y2=verts(face_vert(i,2),2);
                x3=faces(i,1);
                y3=faces(i,2);
                %
                % Pontos de Gauss %
                [x_gauss,ngauss]=gausspoints(x1,y1,x2,y2,x3,y3,'1D',6);
                %
                flux_phi_faces(i,1)=0;
                flux_phi_faces(i,2)=0;
                %
                for j=1:ngauss
                    aux_x=SolutionDiffusion(solution,x_gauss(j,1),x_gauss(j,2),'xflux')*x_gauss(j,3);
                    aux_y=SolutionDiffusion(solution,x_gauss(j,1),x_gauss(j,2),'yflux')*x_gauss(j,3);
                    %
                    flux_phi_faces(i,1)=flux_phi_faces(i,1)+aux_x;
                    flux_phi_faces(i,2)=flux_phi_faces(i,2)+aux_y;
                end
            elseif metodo=='WLS_8'
                %
                % 8ª Ordem %
                x1=verts(face_vert(i,1),1);
                y1=verts(face_vert(i,1),2);
                x2=verts(face_vert(i,2),1);
                y2=verts(face_vert(i,2),2);
                x3=faces(i,1);
                y3=faces(i,2);
                %
                % Pontos de Gauss %
                [x_gauss,ngauss]=gausspoints(x1,y1,x2,y2,x3,y3,'1D',8);
                %
                flux_phi_faces(i,1)=0;
                flux_phi_faces(i,2)=0;
                %
                for j=1:ngauss
                    aux_x=SolutionDiffusion(solution,x_gauss(j,1),x_gauss(j,2),'xflux');
                    aux_y=SolutionDiffusion(solution,x_gauss(j,1),x_gauss(j,2),'yflux');
                    flux_phi_faces(i,1)=flux_phi_faces(i,1)+aux_x*face_area(i,2);
                    flux_phi_faces(i,2)=flux_phi_faces(i,2)+aux_y*face_area(i,3);
                end
            else
                error('\n\nERRO:Não Implementado\n\n');
            end
        else
            error('\n\nERRO:Não Implementado\n\n');
        end
    end
end
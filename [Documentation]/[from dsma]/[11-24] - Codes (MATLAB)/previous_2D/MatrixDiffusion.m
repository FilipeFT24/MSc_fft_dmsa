function [A,source,source_faces,source_cells]=MatrixDiffusion(dirichlet,order,T, T_border,  constrained_source,G,stencil_cells,stencil_faces,stencil_size)
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
global vert_cells vert_cell_num vert_face_num u_convec_x u_convec_y gamma_diff;
global phi lap_phi  phi_faces flux_phi_faces face_w_gauss solution restos mix_method neuman increase_gauss_points dimensional_correction robin;
%
% Inicialização da Matriz A -> Não esquecer que no fim tenho de converter a Matriz A para uma matriz
% esparsa %
% A=zeros(cell_num);
A=sparse(cell_num,cell_num);
contador_aux=0;
for i=1:cell_num
    if rem(i*10,cell_num)==0
        fprintf(' %d%% ',i*100/cell_num);
    end
    %
    % Celula Selecionada %
    cell=i;
    source_faces(cell)=0;
    %
    if dirichlet
        %
        % Faces da Célula %
        for j=1:cell_face_num(cell)
            %
            
            
            % Face Selecionada %
            face=cell_faces(cell,j);
            
            check2=0;
            check4=0;
            check6=0;
            check8=0;
            check = 0;
            
            % Coordenadas da Face
            xf=faces(face,1);
            yf=faces(face,2);
            
            % Area da Face %
            area=face_area(face,1);
            
            % Normal Exterior da Face %
            nx=cell_norm(cell,2*j-1);
            ny=cell_norm(cell,2*j);
            
            % Pontos de Gauss da Face %
%             if order==2
                x1=xf;
                y1=yf;
                P1=1;
%             elseif order==4
%                 
%                 if increase_gauss_points
%                     % Ponto de Gauss 1 %
%                     x1=G{face}(1,1);
%                     y1=G{face}(1,2);
%                     P1=G{face}(1,3);
%                     %
%                     % Ponto de Gauss 2 %
%                     x2=G{face}(2,1);
%                     y2=G{face}(2,2);
%                     P2=G{face}(2,3);
%                     %
%                     % Ponto de Gauss 3 %
%                     x3=G{face}(3,1);
%                     y3=G{face}(3,2);
%                     P3=G{face}(3,3);
%                     
%                 else
%                     % Ponto de Gauss 1 %
%                     x1=G{face}(1,1);
%                     y1=G{face}(1,2);
%                     P1=G{face}(1,3);
%                     %
%                     % Ponto de Gauss 2 %
%                     x2=G{face}(2,1);
%                     y2=G{face}(2,2);
%                     P2=G{face}(2,3);
%                     %
%                 end
%             elseif order==6
%                 % Ponto de Gauss 1 %
%                 x1=G{face}(1,1);
%                 y1=G{face}(1,2);
%                 P1=G{face}(1,3);
%                 %
%                 % Ponto de Gauss 2 %
%                 x2=G{face}(2,1);
%                 y2=G{face}(2,2);
%                 P2=G{face}(2,3);
%                 %
%                 % Ponto de Gauss 3 %
%                 x3=G{face}(3,1);
%                 y3=G{face}(3,2);
%                 P3=G{face}(3,3);
%             elseif order==8
%                 % Ponto de Gauss 1 %
%                 x1=G{face}(1,1);
%                 y1=G{face}(1,2);
%                 P1=G{face}(1,3);
%                 %
%                 % Ponto de Gauss 2 %
%                 x2=G{face}(2,1);
%                 y2=G{face}(2,2);
%                 P2=G{face}(2,3);
%                 %
%                 % Ponto de Gauss 3 %
%                 x3=G{face}(3,1);
%                 y3=G{face}(3,2);
%                 P3=G{face}(3,3);
%                 %
%                 % Ponto de Gauss 4 %
%                 x4=G{face}(4,1);
%                 y4=G{face}(4,2);
%                 P4=G{face}(4,3);
%             end
            %
            for k=1:stencil_size(face,2)
                %
                % Celulas do Stencil %
                cell_aux=stencil_cells(face,k);
                %
                % Constantes %
                C=T{face}(:,k);
                %
                %                 if order==2
                %
                % 2ª Ordem %
                
                if face_w_gauss && face_bound(face,1)==1 && ~mix_method && ~increase_gauss_points && (face_bound(face,4)==0 || face_bound(face,4)==1 || face_bound(face,4)==2)
                    C=T_border{face}(:,k);
                    aux_x1=PolyReconstruction2ndOrder(C,x1,y1,xf,yf,'xflux',face);
                    aux_y1=PolyReconstruction2ndOrder(C,x1,y1,xf,yf,'yflux',face);
                else
                    aux_x1=PolyReconstruction2ndOrder(C,x1,y1,xf,yf,'xflux',face);
                    aux_y1=PolyReconstruction2ndOrder(C,x1,y1,xf,yf,'yflux',face);
                end
                %
                aux_x=aux_x1*area*nx;
                aux_y=aux_y1*area*ny;
                %                 elseif order==4
                %                     %
                %                     % 4ª Ordem %
                %
                %                     if face_w_gauss && face_bound(face,1)==1 && ~mix_method && ~increase_gauss_points && (face_bound(face,4)==0 || face_bound(face,4)==1 || face_bound(face,4)==2)
                %
                %                     C=T_border{face}(:,k);
                %                     aux_x1=PolyReconstruction4thOrder(C,x1,y1,xf,yf,'xflux',face);
                %
                %                     aux_x2=PolyReconstruction4thOrder(C,x2,y2,xf,yf,'xflux',face);
                %
                %                     %
                %
                %                     aux_y1=PolyReconstruction4thOrder(C,x1,y1,xf,yf,'yflux',face);
                %
                %                     aux_y2=PolyReconstruction4thOrder(C,x2,y2,xf,yf,'yflux',face);
                %
                %                     aux_x=(aux_x1*P1+aux_x2*P2)*area*nx;
                %                     aux_y=(aux_y1*P1+aux_y2*P2)*area*ny;
                %
                %                     elseif increase_gauss_points
                %
                %                     aux_x1=PolyReconstruction4thOrder(C,x1,y1,xf,yf,'xflux',face);
                %                     aux_x2=PolyReconstruction4thOrder(C,x2,y2,xf,yf,'xflux',face);
                %                     aux_x3=PolyReconstruction4thOrder(C,x3,y3,xf,yf,'xflux',face);
                %                     %
                %                     aux_y1=PolyReconstruction4thOrder(C,x1,y1,xf,yf,'yflux',face);
                %                     aux_y2=PolyReconstruction4thOrder(C,x2,y2,xf,yf,'yflux',face);
                %                     aux_y3=PolyReconstruction4thOrder(C,x3,y3,xf,yf,'yflux',face);
                %
                %                     %
                %                     aux_x=(aux_x1*P1+aux_x2*P2+aux_x3*P3)*area*nx;
                %                     aux_y=(aux_y1*P1+aux_y2*P2+aux_y3*P3)*area*ny;
                %
                %
                %                     else%if (neuman && ~face_bound(face,4)) || ~neuman
                %
                %                     aux_x1=PolyReconstruction4thOrder(C,x1,y1,xf,yf,'xflux',face);
                %                     aux_x2=PolyReconstruction4thOrder(C,x2,y2,xf,yf,'xflux',face);
                %                     %
                %                     aux_y1=PolyReconstruction4thOrder(C,x1,y1,xf,yf,'yflux',face);
                %                     aux_y2=PolyReconstruction4thOrder(C,x2,y2,xf,yf,'yflux',face);
                %
                %                     aux_x=(aux_x1*P1+aux_x2*P2)*area*nx;
                %                     aux_y=(aux_y1*P1+aux_y2*P2)*area*ny;
                %
                %
                %                     end
                %
                %                 elseif order==6
                %                     %
                %                     % 6ª Ordem %
                %                     if face_w_gauss && face_bound(face,1)==1  && ~mix_method && ~increase_gauss_points && (face_bound(face,4)==0 || face_bound(face,4)==1 || face_bound(face,4)==2)
                %                     C=T_border{face}(:,k);
                %                     aux_x1=PolyReconstruction6thOrder(C,x1,y1,xf,yf,'xflux',face);
                %
                %                     aux_x2=PolyReconstruction6thOrder(C,x2,y2,xf,yf,'xflux',face);
                %
                %                     aux_x3=PolyReconstruction6thOrder(C,x3,y3,xf,yf,'xflux',face);
                %                     %
                %
                %                     aux_y1=PolyReconstruction6thOrder(C,x1,y1,xf,yf,'yflux',face);
                %
                %                     aux_y2=PolyReconstruction6thOrder(C,x2,y2,xf,yf,'yflux',face);
                %
                %                     aux_y3=PolyReconstruction6thOrder(C,x3,y3,xf,yf,'yflux',face);
                %                     else
                %                     aux_x1=PolyReconstruction6thOrder(C,x1,y1,xf,yf,'xflux',face);
                %                     aux_x2=PolyReconstruction6thOrder(C,x2,y2,xf,yf,'xflux',face);
                %                     aux_x3=PolyReconstruction6thOrder(C,x3,y3,xf,yf,'xflux',face);
                %                     %
                %                     aux_y1=PolyReconstruction6thOrder(C,x1,y1,xf,yf,'yflux',face);
                %                     aux_y2=PolyReconstruction6thOrder(C,x2,y2,xf,yf,'yflux',face);
                %                     aux_y3=PolyReconstruction6thOrder(C,x3,y3,xf,yf,'yflux',face);
                %                      end
                %                     %
                %                     aux_x=(aux_x1*P1+aux_x2*P2+aux_x3*P3)*area*nx;
                %                     aux_y=(aux_y1*P1+aux_y2*P2+aux_y3*P3)*area*ny;
                %
                %                 elseif order==8
                %                     %
                %                     % 8ª Ordem
                %                     if face_w_gauss && face_bound(face,1)==1 && ~mix_method && ~increase_gauss_points && (face_bound(face,4)==0 || face_bound(face,4)==1 || face_bound(face,4)==2)
                %                     C=T_border{face}(:,k);
                %                     aux_x1=PolyReconstruction8thOrder(C,x1,y1,xf,yf,'xflux',face);
                %                     aux_x2=PolyReconstruction8thOrder(C,x2,y2,xf,yf,'xflux',face);
                %                     aux_x3=PolyReconstruction8thOrder(C,x3,y3,xf,yf,'xflux',face);
                %                     aux_x4=PolyReconstruction8thOrder(C,x4,y4,xf,yf,'xflux',face);
                %                     %
                %                     aux_y1=PolyReconstruction8thOrder(C,x1,y1,xf,yf,'yflux',face);
                %                     aux_y2=PolyReconstruction8thOrder(C,x2,y2,xf,yf,'yflux',face);
                %                     aux_y3=PolyReconstruction8thOrder(C,x3,y3,xf,yf,'yflux',face);
                %                     aux_y4=PolyReconstruction8thOrder(C,x4,y4,xf,yf,'yflux',face);
                %
                %                     else
                %                     aux_x1=PolyReconstruction8thOrder(C,x1,y1,xf,yf,'xflux',face);
                %                     aux_x2=PolyReconstruction8thOrder(C,x2,y2,xf,yf,'xflux',face);
                %                     aux_x3=PolyReconstruction8thOrder(C,x3,y3,xf,yf,'xflux',face);
                %                     aux_x4=PolyReconstruction8thOrder(C,x4,y4,xf,yf,'xflux',face);
                %                     %
                %                     aux_y1=PolyReconstruction8thOrder(C,x1,y1,xf,yf,'yflux',face);
                %                     aux_y2=PolyReconstruction8thOrder(C,x2,y2,xf,yf,'yflux',face);
                %                     aux_y3=PolyReconstruction8thOrder(C,x3,y3,xf,yf,'yflux',face);
                %                     aux_y4=PolyReconstruction8thOrder(C,x4,y4,xf,yf,'yflux',face);
                %                     end
                %                     aux_x=(aux_x1*P1+aux_x2*P2+aux_x3*P3+aux_x4*P4)*area*nx;
                %                     aux_y=(aux_y1*P1+aux_y2*P2+aux_y3*P3+aux_y4*P4)*area*ny;
                %                 else
                %                     error('\n\nERRO: Ordem Não Implementada\n\n');
                %                 end
                %
                A(cell,cell_aux)=A(cell,cell_aux)+(aux_x+aux_y);
            end
            %
            z=stencil_size(face,2);
            
            
            
            ajuda1=-1;
            factor=0;
            if face_w_gauss && face_bound(face,1)==1 && ~mix_method && (face_bound(face,4)==1 || face_bound(face,4)==0 || face_bound(face,4)==2)
                for ko=1:stencil_size(face,3)
                    if stencil_faces(face,ko)==face
                        ajuda1=ko;
                    end
                end
                
                if ajuda1 ~=-1
                    
                    stencil_faces2=stencil_faces(face,:);
                    stencil_faces2(ajuda1)=[];
                    
                    stencil_faces(face,:)=[stencil_faces2 0];
                    factor=1;
                end
            end
            
            for k=1:stencil_size(face,3)-1*factor
                
                % Face do Stencil %
                face_aux=stencil_faces(face,k);
                
                % Constantes %
                C=T{face}(:,z+k);
                
%                 if order==2
                    %
                    if face_w_gauss && face_bound(face,1)==1 && ~mix_method && face_bound(face_aux,4)==0
                        C=T_border{face}(:,z+k);
                        aux_x1=PolyReconstruction2ndOrder(C,x1,y1,xf,yf,'xflux',face);
                        aux_y1=PolyReconstruction2ndOrder(C,x1,y1,xf,yf,'yflux',face);
                        
                        
                        if face_aux ~=face
                            aux_x=(aux_x1*P1)*area*nx*phi_faces(face_aux);
                            aux_y=(aux_y1*P1)*area*ny*phi_faces(face_aux);
                            
                            if check2==0  && face_w_gauss && face_bound(face,4)==0
                                aux_1xs=PolyReconstruction2ndOrder(constrained_source{face},x1,y1,xf,yf,'xflux',face);
                                aux_1ys=PolyReconstruction2ndOrder(constrained_source{face},x1,y1,xf,yf,'yflux',face);
                                
                                aux_x=aux_x+(aux_1xs*P1)*area*nx;
                                aux_y=aux_y+(aux_1ys*P1)*area*ny;
                                check2=1;
                            end
                        end
                        
                    else
                        if face_w_gauss && face_bound(face,1)==1 && (face_bound(face_aux,4)==1 || face_bound(face_aux,4)==2)
                            C=T_border{face}(:,z+k);
                        end
                        aux_x1=PolyReconstruction2ndOrder(C,x1,y1,xf,yf,'xflux',face);
                        aux_y1=PolyReconstruction2ndOrder(C,x1,y1,xf,yf,'yflux',face);
                        %

                        if neuman && face_bound(face_aux,4)==1
                            
                            normal_x=1*face_area(face_aux,1);
                            normal_y=1*face_area(face_aux,1);
                            
                            
                            if ~dimensional_correction
                                normal_x=normal_x/face_area(face_aux,1);
                                normal_y=normal_y/face_area(face_aux,1);
                            end
                            
                            
                            if face_aux ~=face || ~face_w_gauss
                                aux_x=aux_x1*area*nx*normal_x*SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux');
                                aux_y=aux_y1*area*ny*normal_y*SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux');
                                if check2==0 && face_w_gauss && face_bound(face,4)==1
                                    aux_1xs=PolyReconstruction2ndOrder(constrained_source{face},x1,y1,xf,yf,'xflux',face);
                                    aux_1ys=PolyReconstruction2ndOrder(constrained_source{face},x1,y1,xf,yf,'yflux',face);
                                    
                                    aux_x=aux_x+(aux_1xs*P1)*area*nx;
                                    aux_y=aux_y+(aux_1ys*P1)*area*ny;
                                    check2=1;
                                end
                            end
                            
                        elseif robin && face_bound(face_aux,4)==2
                            
                            normal_x=1*face_area(face_aux,1);
                            normal_y=1*face_area(face_aux,1);
                            
                            
                            if ~dimensional_correction
                                normal_x=normal_x/face_area(face_aux,1);
                                normal_y=normal_y/face_area(face_aux,1);
                            end
                            
                            if face_aux ~=face || ~face_w_gauss
                                
                                aux_x=aux_x1*area*nx*normal_x*(gamma_diff/(u_convec_x))*SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux');
                                aux_y=aux_y1*area*ny*normal_y*(gamma_diff/(u_convec_x))*SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux');
                                aux_x=aux_x+aux_x1*area*normal_x*nx*phi_faces(face_aux);
                                aux_y=aux_y+aux_y1*area*normal_y*ny*phi_faces(face_aux);
                                
                                if check2==0  && face_w_gauss && face_bound(face,4)==2
                                    aux_1xs=PolyReconstruction2ndOrder(constrained_source{face},x1,y1,xf,yf,'xflux',face);
                                    aux_1ys=PolyReconstruction2ndOrder(constrained_source{face},x1,y1,xf,yf,'yflux',face);
                                    
                                    aux_x=aux_x+(aux_1xs*P1)*area*nx;
                                    aux_y=aux_y+(aux_1ys*P1)*area*ny;
                                    check2=1;
                                end
                            end
                            
                        else
                            aux_x=aux_x1*area*nx*phi_faces(face_aux);
                            aux_y=aux_y1*area*ny*phi_faces(face_aux);
                        end
                    end
                    
%                 elseif order==4
%                     %
%                     % 4ª Ordem %
%                     
%                     if face_w_gauss && face_bound(face,1)==1 && ~mix_method && ~increase_gauss_points && face_bound(face_aux,4)==0
%                         
%                         C=T_border{face}(:,z+k);
%                         aux_x1=PolyReconstruction4thOrder(C,x1,y1,xf,yf,'xflux',face);
%                         
%                         aux_x2=PolyReconstruction4thOrder(C,x2,y2,xf,yf,'xflux',face);
%                         
%                         %
%                         
%                         aux_y1=PolyReconstruction4thOrder(C,x1,y1,xf,yf,'yflux',face);
%                         
%                         aux_y2=PolyReconstruction4thOrder(C,x2,y2,xf,yf,'yflux',face);
%                         
%                         
%                         if face_aux ~=face
%                             aux_x=(aux_x1*P1+aux_x2*P2)*area*nx*phi_faces(face_aux);
%                             aux_y=(aux_y1*P1+aux_y2*P2)*area*ny*phi_faces(face_aux);
%                             
%                             if check4==0  && face_w_gauss && face_bound(face,4)==0
%                                 aux_1xs=PolyReconstruction4thOrder(constrained_source{face},x1,y1,xf,yf,'xflux',face);
%                                 aux_2xs=PolyReconstruction4thOrder(constrained_source{face},x2,y2,xf,yf,'xflux',face);
%                                 
%                                 
%                                 aux_1ys=PolyReconstruction4thOrder(constrained_source{face},x1,y1,xf,yf,'yflux',face);
%                                 aux_2ys=PolyReconstruction4thOrder(constrained_source{face},x2,y2,xf,yf,'yflux',face);
%                                 
%                                 
%                                 aux_x=aux_x+(aux_1xs*P1+aux_2xs*P2)*area*nx;
%                                 aux_y=aux_y+(aux_1ys*P1+aux_2ys*P2)*area*ny;
%                                 check4=1;
%                             end
%                         end
%                         
%                     elseif increase_gauss_points
%                         aux_x1=PolyReconstruction4thOrder(C,x1,y1,xf,yf,'xflux',face);
%                         aux_x2=PolyReconstruction4thOrder(C,x2,y2,xf,yf,'xflux',face);
%                         aux_x3=PolyReconstruction4thOrder(C,x3,y3,xf,yf,'xflux',face);
%                         %
%                         aux_y1=PolyReconstruction4thOrder(C,x1,y1,xf,yf,'yflux',face);
%                         aux_y2=PolyReconstruction4thOrder(C,x2,y2,xf,yf,'yflux',face);
%                         aux_y3=PolyReconstruction4thOrder(C,x3,y3,xf,yf,'yflux',face);
%                         
%                         %
%                         aux_x=(aux_x1*P1+aux_x2*P2+aux_x3*P3)*area*nx*phi_faces(face_aux);
%                         aux_y=(aux_y1*P1+aux_y2*P2+aux_y3*P3)*area*ny*phi_faces(face_aux);
%                         
%                         
%                     else
%                         if face_w_gauss && face_bound(face,1)==1 && (face_bound(face_aux,4)==1 || face_bound(face_aux,4)==2)
%                             C=T_border{face}(:,z+k);
%                         end
%                         aux_x1=PolyReconstruction4thOrder(C,x1,y1,xf,yf,'xflux',face);
%                         aux_x2=PolyReconstruction4thOrder(C,x2,y2,xf,yf,'xflux',face);
%                         %
%                         aux_y1=PolyReconstruction4thOrder(C,x1,y1,xf,yf,'yflux',face);
%                         aux_y2=PolyReconstruction4thOrder(C,x2,y2,xf,yf,'yflux',face);
%                         %
%                         if neuman && face_bound(face_aux,4)==1
%                             
%                             normal_x=1*face_area(face_aux,1);
%                             normal_y=1*face_area(face_aux,1);
%                             
%                             if ~dimensional_correction
%                                 normal_x=normal_x/face_area(face_aux,1);
%                                 normal_y=normal_y/face_area(face_aux,1);
%                             end
%                             
%                             if face_aux ~=face || ~face_w_gauss
%                                 
%                                 aux_x=(aux_x1*P1+aux_x2*P2)*area*nx*normal_x*(SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux'));
%                                 aux_y=(aux_y1*P1+aux_y2*P2)*area*ny*normal_y*(SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux'));
%                                 
%                                 if check4==0 && face_w_gauss && face_bound(face,4)==1
%                                     contador_aux=contador_aux+1;
%                                     aux_1xs=PolyReconstruction4thOrder(constrained_source{face},x1,y1,xf,yf,'xflux',face);
%                                     aux_2xs=PolyReconstruction4thOrder(constrained_source{face},x2,y2,xf,yf,'xflux',face);
%                                     
%                                     
%                                     aux_1ys=PolyReconstruction4thOrder(constrained_source{face},x1,y1,xf,yf,'yflux',face);
%                                     aux_2ys=PolyReconstruction4thOrder(constrained_source{face},x2,y2,xf,yf,'yflux',face);
%                                     
%                                     
%                                     aux_x=aux_x+(aux_1xs*P1+aux_2xs*P2)*area*nx;
%                                     aux_y=aux_y+(aux_1ys*P1+aux_2ys*P2)*area*ny;
%                                     check4=1;
%                                 end
%                             end
%                             
%                         elseif robin && face_bound(face_aux,4)==2
%                             
%                             normal_x=1*face_area(face_aux,1);
%                             normal_y=1*face_area(face_aux,1);
%                             
%                             
%                             if ~dimensional_correction
%                                 normal_x=normal_x/face_area(face_aux,1);
%                                 normal_y=normal_y/face_area(face_aux,1);
%                             end
%                             
%                             if face_aux ~=face || ~face_w_gauss
%                                 aux_x=(aux_x1*P1+aux_x2*P2)*area*nx*normal_x*(gamma_diff/(u_convec_x))*(SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux'));
%                                 aux_y=(aux_y1*P1+aux_y2*P2)*area*ny*normal_y*(gamma_diff/(u_convec_x))*(SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux'));
%                                 
%                                 aux_x=aux_x+(aux_x1*P1+aux_x2*P2)*area*normal_x*nx*phi_faces(face_aux);
%                                 aux_y=aux_y+(aux_y1*P1+aux_y2*P2)*area*normal_y*ny*phi_faces(face_aux);
%                                 
%                                 if check4==0 && face_w_gauss && face_bound(face,4)==2
%                                     contador_aux=contador_aux+1;
%                                     aux_1xs=PolyReconstruction4thOrder(constrained_source{face},x1,y1,xf,yf,'xflux',face);
%                                     aux_2xs=PolyReconstruction4thOrder(constrained_source{face},x2,y2,xf,yf,'xflux',face);
%                                     
%                                     
%                                     aux_1ys=PolyReconstruction4thOrder(constrained_source{face},x1,y1,xf,yf,'yflux',face);
%                                     aux_2ys=PolyReconstruction4thOrder(constrained_source{face},x2,y2,xf,yf,'yflux',face);
%                                     
%                                     
%                                     aux_x=aux_x+(aux_1xs*P1+aux_2xs*P2)*area*nx;
%                                     aux_y=aux_y+(aux_1ys*P1+aux_2ys*P2)*area*ny;
%                                     check4=1;
%                                 end
%                             end
%                             
%                         else
%                             aux_x=(aux_x1*P1+aux_x2*P2)*area*nx*phi_faces(face_aux);
%                             aux_y=(aux_y1*P1+aux_y2*P2)*area*ny*phi_faces(face_aux);
%                         end
%                         
%                         
%                         if abs(aux_x)>1 || abs(aux_y)>1
%                             pause(0.1)
%                         end
%                         
%                         
%                     end
%                 elseif order==6
%                     %
%                     % 6ª Ordem %
%                     if face_w_gauss && face_bound(face,1)==1 && ~mix_method && face_bound(face_aux,4)==0
%                         C=T_border{face}(:,z+k);
%                         aux_x1=PolyReconstruction6thOrder(C,x1,y1,xf,yf,'xflux',face);
%                         
%                         aux_x2=PolyReconstruction6thOrder(C,x2,y2,xf,yf,'xflux',face);
%                         
%                         aux_x3=PolyReconstruction6thOrder(C,x3,y3,xf,yf,'xflux',face);
%                         
%                         
%                         aux_y1=PolyReconstruction6thOrder(C,x1,y1,xf,yf,'yflux',face);
%                         
%                         aux_y2=PolyReconstruction6thOrder(C,x2,y2,xf,yf,'yflux',face);
%                         
%                         aux_y3=PolyReconstruction6thOrder(C,x3,y3,xf,yf,'yflux',face);
%                         
%                         
%                         if face_aux ~=face
%                             aux_x=(aux_x1*P1+aux_x2*P2+aux_x3*P3)*area*nx*phi_faces(face_aux);
%                             aux_y=(aux_y1*P1+aux_y2*P2+aux_y3*P3)*area*ny*phi_faces(face_aux);
%                             
%                             if check6==0  && face_w_gauss && face_bound(face,4)==0
%                                 aux_1xs=PolyReconstruction6thOrder(constrained_source{face},x1,y1,xf,yf,'xflux',face);
%                                 aux_2xs=PolyReconstruction6thOrder(constrained_source{face},x2,y2,xf,yf,'xflux',face);
%                                 aux_3xs=PolyReconstruction6thOrder(constrained_source{face},x3,y3,xf,yf,'xflux',face);
%                                 
%                                 aux_1ys=PolyReconstruction6thOrder(constrained_source{face},x1,y1,xf,yf,'yflux',face);
%                                 aux_2ys=PolyReconstruction6thOrder(constrained_source{face},x2,y2,xf,yf,'yflux',face);
%                                 aux_3ys=PolyReconstruction6thOrder(constrained_source{face},x3,y3,xf,yf,'yflux',face);
%                                 
%                                 aux_x=aux_x+(aux_1xs*P1+aux_2xs*P2+aux_3xs*P3)*area*nx;
%                                 aux_y=aux_y+(aux_1ys*P1+aux_2ys*P2+aux_3ys*P3)*area*ny;
%                                 check6=1;
%                             end
%                         end
%                         
%                     else
%                         if face_w_gauss && face_bound(face,1)==1 && (face_bound(face_aux,4)==1 || face_bound(face_aux,4)==2)
%                             C=T_border{face}(:,z+k);
%                         end
%                         aux_x1=PolyReconstruction6thOrder(C,x1,y1,xf,yf,'xflux',face);
%                         aux_x2=PolyReconstruction6thOrder(C,x2,y2,xf,yf,'xflux',face);
%                         aux_x3=PolyReconstruction6thOrder(C,x3,y3,xf,yf,'xflux',face);
%                         %
%                         aux_y1=PolyReconstruction6thOrder(C,x1,y1,xf,yf,'yflux',face);
%                         aux_y2=PolyReconstruction6thOrder(C,x2,y2,xf,yf,'yflux',face);
%                         aux_y3=PolyReconstruction6thOrder(C,x3,y3,xf,yf,'yflux',face);
%                         
%                         if neuman && face_bound(face_aux,4)==1
%                             normal_x=1*face_area(face_aux,1);
%                             normal_y=1*face_area(face_aux,1);
%                             
%                             if ~dimensional_correction
%                                 normal_x=normal_x/face_area(face_aux,1);
%                                 normal_y=normal_y/face_area(face_aux,1);
%                             end
%                             
%                             if face_aux ~=face || ~face_w_gauss
%                                 aux_x=(aux_x1*P1+aux_x2*P2+aux_x3*P3)*area*nx*normal_x*(SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux'));
%                                 aux_y=(aux_y1*P1+aux_y2*P2+aux_y3*P3)*area*ny*normal_y*(SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux'));
%                                 if check6==0 && face_w_gauss && face_bound(face,4)==1
%                                     aux_1xs=PolyReconstruction6thOrder(constrained_source{face},x1,y1,xf,yf,'xflux',face);
%                                     aux_2xs=PolyReconstruction6thOrder(constrained_source{face},x2,y2,xf,yf,'xflux',face);
%                                     aux_3xs=PolyReconstruction6thOrder(constrained_source{face},x3,y3,xf,yf,'xflux',face);
%                                     
%                                     aux_1ys=PolyReconstruction6thOrder(constrained_source{face},x1,y1,xf,yf,'yflux',face);
%                                     aux_2ys=PolyReconstruction6thOrder(constrained_source{face},x2,y2,xf,yf,'yflux',face);
%                                     aux_3ys=PolyReconstruction6thOrder(constrained_source{face},x3,y3,xf,yf,'yflux',face);
%                                     
%                                     aux_x=aux_x+(aux_1xs*P1+aux_2xs*P2+aux_3xs*P3)*area*nx;
%                                     aux_y=aux_y+(aux_1ys*P1+aux_2ys*P2+aux_3ys*P3)*area*ny;
%                                     check6=1;
%                                 end
%                             end
%                             
%                             
%                         elseif robin && face_bound(face_aux,4)==2
%                             normal_x=1*face_area(face_aux,1);
%                             normal_y=1*face_area(face_aux,1);
%                             
%                             
%                             if ~dimensional_correction
%                                 normal_x=normal_x/face_area(face_aux,1);
%                                 normal_y=normal_y/face_area(face_aux,1);
%                             end
%                             
%                             if face_aux ~=face || ~face_w_gauss
%                                 aux_x=(aux_x1*P1+aux_x2*P2+aux_x3*P3)*area*nx*normal_x*(gamma_diff/(u_convec_x))*(SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux'));
%                                 aux_y=(aux_y1*P1+aux_y2*P2+aux_y3*P3)*area*ny*normal_y*(gamma_diff/(u_convec_x))*(SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux'));
%                                 
%                                 aux_x=aux_x+(aux_x1*P1+aux_x2*P2+aux_x3*P3)*normal_x*area*nx*phi_faces(face_aux);
%                                 aux_y=aux_y+(aux_y1*P1+aux_y2*P2+aux_y3*P3)*normal_y*area*ny*phi_faces(face_aux);
%                                 
%                                 if check6==0  && face_w_gauss && face_bound(face,4)==2
%                                     aux_1xs=PolyReconstruction6thOrder(constrained_source{face},x1,y1,xf,yf,'xflux',face);
%                                     aux_2xs=PolyReconstruction6thOrder(constrained_source{face},x2,y2,xf,yf,'xflux',face);
%                                     aux_3xs=PolyReconstruction6thOrder(constrained_source{face},x3,y3,xf,yf,'xflux',face);
%                                     
%                                     aux_1ys=PolyReconstruction6thOrder(constrained_source{face},x1,y1,xf,yf,'yflux',face);
%                                     aux_2ys=PolyReconstruction6thOrder(constrained_source{face},x2,y2,xf,yf,'yflux',face);
%                                     aux_3ys=PolyReconstruction6thOrder(constrained_source{face},x3,y3,xf,yf,'yflux',face);
%                                     
%                                     aux_x=aux_x+(aux_1xs*P1+aux_2xs*P2+aux_3xs*P3)*area*nx;
%                                     aux_y=aux_y+(aux_1ys*P1+aux_2ys*P2+aux_3ys*P3)*area*ny;
%                                     check6=1;
%                                 end
%                             end
%                             
%                         else
%                             aux_x=(aux_x1*P1+aux_x2*P2+aux_x3*P3)*area*nx*phi_faces(face_aux);
%                             aux_y=(aux_y1*P1+aux_y2*P2+aux_y3*P3)*area*ny*phi_faces(face_aux);
%                         end
%                     end
%                     
%                     
%                 elseif order==8
%                     %
%                     % 8ª Ordem
%                     if face_w_gauss && face_bound(face,1)==1 && ~mix_method && face_bound(face_aux,4)==0
%                         C=T_border{face}(:,z+k);
%                         aux_x1=PolyReconstruction8thOrder(C,x1,y1,xf,yf,'xflux',face);
%                         
%                         aux_x2=PolyReconstruction8thOrder(C,x2,y2,xf,yf,'xflux',face);
%                         
%                         aux_x3=PolyReconstruction8thOrder(C,x3,y3,xf,yf,'xflux',face);
%                         aux_x4=PolyReconstruction8thOrder(C,x4,y4,xf,yf,'xflux',face);
%                         
%                         %
%                         
%                         aux_y1=PolyReconstruction8thOrder(C,x1,y1,xf,yf,'yflux',face);
%                         
%                         aux_y2=PolyReconstruction8thOrder(C,x2,y2,xf,yf,'yflux',face);
%                         
%                         aux_y3=PolyReconstruction8thOrder(C,x3,y3,xf,yf,'yflux',face);
%                         aux_y4=PolyReconstruction8thOrder(C,x4,y4,xf,yf,'yflux',face);
%                         
%                         if face_aux ~=face
%                             aux_x=(aux_x1*P1+aux_x2*P2+aux_x3*P3+aux_x4*P4)*area*nx*phi_faces(face_aux);
%                             aux_y=(aux_y1*P1+aux_y2*P2+aux_y3*P3+aux_y4*P4)*area*ny*phi_faces(face_aux);
%                             
%                             if check8==0  && face_w_gauss && face_bound(face,4)==0
%                                 aux_1xs=PolyReconstruction8thOrder(constrained_source{face},x1,y1,xf,yf,'xflux',face);
%                                 aux_2xs=PolyReconstruction8thOrder(constrained_source{face},x2,y2,xf,yf,'xflux',face);
%                                 aux_3xs=PolyReconstruction8thOrder(constrained_source{face},x3,y3,xf,yf,'xflux',face);
%                                 aux_4xs=PolyReconstruction8thOrder(constrained_source{face},x4,y4,xf,yf,'xflux',face);
%                                 
%                                 aux_1ys=PolyReconstruction8thOrder(constrained_source{face},x1,y1,xf,yf,'yflux',face);
%                                 aux_2ys=PolyReconstruction8thOrder(constrained_source{face},x2,y2,xf,yf,'yflux',face);
%                                 aux_3ys=PolyReconstruction8thOrder(constrained_source{face},x3,y3,xf,yf,'yflux',face);
%                                 aux_4ys=PolyReconstruction8thOrder(constrained_source{face},x4,y4,xf,yf,'yflux',face);
%                                 
%                                 aux_x=aux_x+(aux_1xs*P1+aux_2xs*P2+aux_3xs*P3+aux_4xs*P4)*area*nx;
%                                 aux_y=aux_y+(aux_1ys*P1+aux_2ys*P2+aux_3ys*P3+aux_4ys*P4)*area*ny;
%                                 check8=1;
%                             end
%                         end
%                         
%                     else
%                         if face_w_gauss && face_bound(face,1)==1 && (face_bound(face_aux,4)==1 || face_bound(face_aux,4)==2)
%                             C=T_border{face}(:,z+k);
%                         end
%                         aux_x1=PolyReconstruction8thOrder(C,x1,y1,xf,yf,'xflux',face);
%                         aux_x2=PolyReconstruction8thOrder(C,x2,y2,xf,yf,'xflux',face);
%                         aux_x3=PolyReconstruction8thOrder(C,x3,y3,xf,yf,'xflux',face);
%                         aux_x4=PolyReconstruction8thOrder(C,x4,y4,xf,yf,'xflux',face);
%                         %
%                         aux_y1=PolyReconstruction8thOrder(C,x1,y1,xf,yf,'yflux',face);
%                         aux_y2=PolyReconstruction8thOrder(C,x2,y2,xf,yf,'yflux',face);
%                         aux_y3=PolyReconstruction8thOrder(C,x3,y3,xf,yf,'yflux',face);
%                         aux_y4=PolyReconstruction8thOrder(C,x4,y4,xf,yf,'yflux',face);
%                         %
%                         if neuman && face_bound(face_aux,4)==1
%                             normal_x=1*face_area(face_aux,1);
%                             normal_y=1*face_area(face_aux,1);
%                             
%                             
%                             if ~dimensional_correction
%                                 normal_x=normal_x/face_area(face_aux,1);
%                                 normal_y=normal_y/face_area(face_aux,1);
%                             end
%                             
%                             if face_aux ~=face || ~face_w_gauss
%                                 aux_x=(aux_x1*P1+aux_x2*P2+aux_x3*P3+aux_x4*P4)*area*nx*normal_x*(SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux'));
%                                 aux_y=(aux_y1*P1+aux_y2*P2+aux_y3*P3+aux_y4*P4)*area*ny*normal_y*(SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux'));
%                                 if check8==0 && face_w_gauss && face_bound(face,4)==1
%                                     aux_1xs=PolyReconstruction8thOrder(constrained_source{face},x1,y1,xf,yf,'xflux',face);
%                                     aux_2xs=PolyReconstruction8thOrder(constrained_source{face},x2,y2,xf,yf,'xflux',face);
%                                     aux_3xs=PolyReconstruction8thOrder(constrained_source{face},x3,y3,xf,yf,'xflux',face);
%                                     aux_4xs=PolyReconstruction8thOrder(constrained_source{face},x4,y4,xf,yf,'xflux',face);
%                                     
%                                     aux_1ys=PolyReconstruction8thOrder(constrained_source{face},x1,y1,xf,yf,'yflux',face);
%                                     aux_2ys=PolyReconstruction8thOrder(constrained_source{face},x2,y2,xf,yf,'yflux',face);
%                                     aux_3ys=PolyReconstruction8thOrder(constrained_source{face},x3,y3,xf,yf,'yflux',face);
%                                     aux_4ys=PolyReconstruction8thOrder(constrained_source{face},x4,y4,xf,yf,'yflux',face);
%                                     
%                                     aux_x=aux_x+(aux_1xs*P1+aux_2xs*P2+aux_3xs*P3+aux_4xs*P4)*area*nx;
%                                     aux_y=aux_y+(aux_1ys*P1+aux_2ys*P2+aux_3ys*P3+aux_4ys*P4)*area*ny;
%                                     check8=1;
%                                 end
%                             end
%                             
%                         elseif robin && face_bound(face_aux,4)==2
%                             normal_x=1*face_area(face_aux,1);
%                             normal_y=1*face_area(face_aux,1);
%                             
%                             
%                             if ~dimensional_correction
%                                 normal_x=normal_x/face_area(face_aux,1);
%                                 normal_y=normal_y/face_area(face_aux,1);
%                             end
%                             
%                             if face_aux ~=face || ~face_w_gauss
%                                 aux_x=(aux_x1*P1+aux_x2*P2+aux_x3*P3+aux_x4*P4)*area*nx*normal_x*(gamma_diff/(u_convec_x))*(SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux'));
%                                 aux_y=(aux_y1*P1+aux_y2*P2+aux_y3*P3+aux_y4*P4)*area*ny*normal_y*(gamma_diff/(u_convec_x))*(SolutionDiffusion(solution,faces(face_aux,1),faces(face_aux,2),'xflux'));
%                                 
%                                 aux_x=aux_x+(aux_x1*P1+aux_x2*P2+aux_x3*P3+aux_x4*P4)*normal_x*area*nx*phi_faces(face_aux);
%                                 aux_y=aux_y+(aux_y1*P1+aux_y2*P2+aux_y3*P3+aux_y4*P4)*normal_y*area*ny*phi_faces(face_aux);
%                                 
%                                 if check8==0  && face_w_gauss && face_bound(face,4)==2
%                                     aux_1xs=PolyReconstruction8thOrder(constrained_source{face},x1,y1,xf,yf,'xflux',face);
%                                     aux_2xs=PolyReconstruction8thOrder(constrained_source{face},x2,y2,xf,yf,'xflux',face);
%                                     aux_3xs=PolyReconstruction8thOrder(constrained_source{face},x3,y3,xf,yf,'xflux',face);
%                                     aux_4xs=PolyReconstruction8thOrder(constrained_source{face},x4,y4,xf,yf,'xflux',face);
%                                     
%                                     aux_1ys=PolyReconstruction8thOrder(constrained_source{face},x1,y1,xf,yf,'yflux',face);
%                                     aux_2ys=PolyReconstruction8thOrder(constrained_source{face},x2,y2,xf,yf,'yflux',face);
%                                     aux_3ys=PolyReconstruction8thOrder(constrained_source{face},x3,y3,xf,yf,'yflux',face);
%                                     aux_4ys=PolyReconstruction8thOrder(constrained_source{face},x4,y4,xf,yf,'yflux',face);
%                                     
%                                     aux_x=aux_x+(aux_1xs*P1+aux_2xs*P2+aux_3xs*P3+aux_4xs*P4)*area*nx;
%                                     aux_y=aux_y+(aux_1ys*P1+aux_2ys*P2+aux_3ys*P3+aux_4ys*P4)*area*ny;
%                                     check8=1;
%                                 end
%                             end
%                             
%                             
%                         else
%                             aux_x=(aux_x1*P1+aux_x2*P2+aux_x3*P3+aux_x4*P4)*area*nx*phi_faces(face_aux);
%                             aux_y=(aux_y1*P1+aux_y2*P2+aux_y3*P3+aux_y4*P4)*area*ny*phi_faces(face_aux);
%                         end
%                     end
%                 else
%                     error('\n\nERRO: Ordem Não Implementada\n\n');
%                 end
                %
                source_faces(cell)=source_faces(cell)+(aux_x+aux_y);
            end
        end
    else
        error('\n\nERRO:Condição de Fronteira Não Implementada\n\n');
    end
    %
    source_cells(cell)=lap_phi(cell)*cell_vol(cell);
    source(cell)=-source_faces(cell);
end
%
fprintf(' ... Fim\n');

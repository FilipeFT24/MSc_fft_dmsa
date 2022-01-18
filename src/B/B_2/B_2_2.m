classdef B_2_2
    methods (Static)
        %% > Wrap-up B_2_2.
        function [pde] = WrapUp_B_2_2(inp,msh,pde)
            % >> ----------------------------------------------------------
            % >> 1.   Assemble matrices Df, Dwf, Pf and Tf.
            %  > 1.1. Compute df array.
            % >> 2.   Assemble matrices A and B.
            %  > 2.1. Assemble matrix A.
            %  > 2.2. Assemble matrix B.
            % >> ----------------------------------------------------------
            % >> Local variables.
            ng = inp.fr.ng;
            V  = [inp.pr.vx,inp.pr.vy];
            G  = [inp.pr.gx,inp.pr.gy];
            
            % >> 1.
            pde = B_2_2.Assemble_Mat_1(msh,pde,ng);
            % >> 2.
            pde = B_2_2.Assemble_Mat_2(msh,pde,V,G);
        end
        
        %% > 1. -----------------------------------------------------------
        function [pde] = Assemble_Mat_1(msh,pde,ng)
            %% > Matrices Df and Dwf.
            % >> Df = [1*{(x-xf)^0}*{(y-yf)^0},1*{(x-xf)^1}*{(y-yf)^0},1*{(x-xf)^0}*{(y-yf)^1},...] = [1,(x-xf),(y-yf),...].
            for i = 1:msh.f.NF
                %  > Face.
                Face(:,i) = [msh.f.mean(1,i),msh.f.mean(2,i)];
                %  > Stencil point coordinates.
                xy_v{i} = msh.s.xy_v_t{i};
                for j = 1:size(xy_v{i},2)
                    %  > x-xf.
                    XY{i}(1,j) = xy_v{i}(1,j)-Face(1,i);
                    XY{i}(2,j) = xy_v{i}(2,j)-Face(2,i);
                    %  > Deal coefficients.
                    for k = 1:pde.pr.numb
                        Df{i}(j,k) = pde.pr.Coeff_1(k).*(XY{i}(1,j).^pde.pr.Exp_1(1,k)).*(XY{i}(2,j).^pde.pr.Exp_1(2,k));
                    end
                end
            end
            
            % >> Dwf = w*Df.
            %  > Boundary faces' cell index.
            for i = 1:size(msh.bnd.f,2)
                bnd_ff(i) = msh.bnd.f{2,i};
                bnd_fc(i) = msh.bnd.f{3,i};
            end
            for i = 1:msh.f.NF
                for j = 1:size(xy_v{i},2)
                    if ~isequal(xy_v{i}(:,j),Face(:,i))
                        %  > Point.
                        Point{i}(:,j) = xy_v{i}(:,j);
                        %  > Distance.
                        d    {i}(j)   = pdist([Point{i}(:,j),Face(:,i)],'euclidean');
                        w    {i}(j)   = pde.wf.wf_1(d{i}(j),2);
                    else
                        %  > Point.
                        %    Remark: If point=face, use cell centroid instead.
                        Point{i}(:,j) = msh.c.mean(:,bnd_fc(bnd_ff == i));
                        %  > Distance.
                        d    {i}(j)   = pdist([Point{i}(:,j),Face(:,i)],'euclidean');
                        w    {i}(j)   = pde.wf.wf_2(d{i}(j),2);
                    end
                    Dwf{i}(j,:) = Df{i}(j,:).*w{i}(j);
                end
            end
            
            %% > Matrices Pf and Tf.
            % >> Pf.
            %  > Pf = inv(Dwf_T*Df)*Dwf_T.
            for i = 1:msh.f.NF
                Inv{i} = eMatrices(transpose(Dwf{i})*Df{i});
                Pf {i} = Inv{i}*transpose(Dwf{i});
            end
            % >> Tf     = ?
            %  > C(1)   = P(1,1)*Phi(1)+P(1,2)*Phi(2)+P(1,3)*Phi(3)+P(1,4)*Phi(4)+...
            %  > C(2)   = P(2,1)*Phi(1)+P(2,2)*Phi(2)+P(2,3)*Phi(3)+P(2,4)*Phi(4)+...
            %  > C(3)   = P(3,1)*Phi(1)+P(3,2)*Phi(2)+P(3,3)*Phi(3)+P(3,4)*Phi(4)+...
            %  > Phi_f  = C(1)+C(2)*(x-xf)+C(3)*(y-yf)+...
            %  > Phi_f  = {P(1,1)*Phi(1)+P(1,2)*Phi(2)+P(1,3)*Phi(3)+P(1,4)*Phi(4)}+
            %           + {P(2,1)*Phi(1)+P(2,2)*Phi(2)+P(2,3)*Phi(3)+P(2,4)*Phi(4)}*(x-xf)+
            %           + {P(3,1)*Phi(1)+P(3,2)*Phi(2)+P(3,3)*Phi(3)+P(3,4)*Phi(4)}*(y-yf)+...
            %  > Phi_f  = {P(1,1)+P(2,1)*(x-xf)+P(3,1)*(y-yf)+...}*Phi(1)+
            %           + {P(1,2)+P(2,2)*(x-xf)+P(3,2)*(y-yf)+...}*Phi(2)+
            %           + {P(1,3)+P(2,3)*(x-xf)+P(3,3)*(y-yf)+...}*Phi(3)+...
            %           + {P(1,4)+P(2,4)*(x-xf)+P(3,4)*(y-yf)+...}*Phi(4)+...
            %
            %  > Tf    -> Equivalent to: [1,(x-xf),(y-yf),...]*Pf*Phi = df*Pf*Phi.
            %    e.g.:    p=1 w/ ns = 8 -> (1x3)*(3x8)*(8x1) = (1x8)*(8x1) = (1x1).
            
            %  > Coordinate transformation function handle.
            [xy_fg,j_fg] = B_1_2.CD_1D();
            
            for i = 1:msh.f.NF
                %  > Face quadrature.
                pde.f.gq{i} = B_1_2.GaussFace_Points(xy_fg,j_fg,ng,msh.f.xy_v{i});
                %  > Convection.
                df      {i} = B_2_2.Compute_df(1,Face(:,i),pde.f.gq{i},pde.pr.Coeff_1,pde.pr.Exp_1);
                Tf_C    {i} = df{i}*Pf{i};
                %  > Diffusion.
                grad_df {i} = B_2_2.Compute_df(2,Face(:,i),pde.f.gq{i},pde.pr.Coeff_2,pde.pr.Exp_2);
                Tf_D    {i} = grad_df{i}*Pf{i};
            end

            % >> Deal...
            %  > Df.
            pde.f.Df   = Df;
            %  > Dwf.
            pde.f.Dwf  = Dwf;
            %  > Pf.
            pde.f.Pf   = Pf;
            %  > Tf.
            pde.f.Tf_C = Tf_C;
            pde.f.Tf_D = Tf_D;
        end
        % >> 1.1. ---------------------------------------------------------
        function [df_ij] = Compute_df(iD,mean_f,fg,Coeff,Exp)
            if iD == 1
                % >> Phi_f.
                for i = 1:size(fg.Points,1)
                    for j = 1:length(Coeff)
                        df_ij(i,j) = Coeff(j).*fg.Weights(i).*((mean_f(1)-fg.Points(i,1)).^Exp(1,j)).*((mean_f(2)-fg.Points(i,2)).^Exp(2,j));
                    end
                end
                df_ij = sum(df_ij,1);
            elseif iD == 2
                % >> grad(Phi_f).
                for i = 1:size(Coeff,1)
                    %  > i=1: X.
                    %  > i=2: Y.
                    for j = 1:size(fg.Points,1)
                        for k = 1:size(Coeff,2)
                            df_ij{i}(j,k) = Coeff(i,k).*fg.Weights(j).*((mean_f(1)-fg.Points(j,1)).^Exp{i}(1,k)).*((mean_f(2)-fg.Points(j,2)).^Exp{i}(2,k));
                        end
                    end
                    df_ij{i} = sum(df_ij{i},1);
                end
                %  > Reshape cell array.
                df_ij = cell2mat(reshape(df_ij,[2,1]));
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [pde] = Assemble_Mat_2(msh,pde,V,G)
            % >> Tf.
            %  > V*Tf_C.
            VxTf_C = cellfun(@(x) x*V(1),pde.f.Tf_C,'un',0);
            VyTf_C = cellfun(@(x) x*V(2),pde.f.Tf_C,'un',0);
            %  > G*Tf_D.
            GxTf_D = cellfun(@(x) x*G(1),pde.f.Tf_D,'un',0);
            GyTf_D = cellfun(@(x) x*G(2),pde.f.Tf_D,'un',0);
            
            % >> Phi_s.
            for i = 1:msh.f.NF
                Phi_s{i} = A_3_2_1.Deal_StencilElem(msh.s.c(:,i));
            end
            
       
            
            
            
            
        end
        % >> 2.2. ---------------------------------------------------------
    end
end
classdef B_2_2
    methods (Static)
        %% > B_2_2.
        function [pde] = WrapUp_B_2_2(inp,msh,pde)
            % >> ----------------------------------------------------------
            % >> 1.   Assemble matrices Df, Dwf, Pf and Tf.
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
        end
        
        %% > 1. -----------------------------------------------------------
        function [pde] = Assemble_Mat_1(msh,pde,ng)
            %% > Matrices Df and Dwf.
            % >> Df.
            %  > Df = [1*{(x-xf)^0}*{(y-yf)^0},1*{(x-xf)^1}*{(y-yf)^0},1*{(x-xf)^0}*{(y-yf)^1},...] = [1,(x-xf),(y-yf),...].
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
            
            % >> Dwf.
            %  > Dwf = w*Df.
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
                        w    {i}(j)   = pde.wf.wf_1(d{i}(j),4);
                    else
                        %  > Point.
                        %    Remark: If point=face, use cell centroid instead.
                        Point{i}(:,j) = msh.c.mean(:,bnd_fc(bnd_ff == i));
                        %  > Distance.
                        d    {i}(j)   = pdist([Point{i}(:,j),Face(:,i)],'euclidean');
                        w    {i}(j)   = pde.wf.wf_2(d{i}(j),4);
                    end
                    Dwf{i}(j,:) = Df{i}(j,:).*w{i}(j);
                end
            end
            
            %% > Matrices Pf and Tf.
            % >> Pf.
            %  > Pf = inv(Dwf_T*Df)*Dwf_T.
            for i = 1:msh.f.NF
                Inv           {i} = eMatrices(transpose(Dwf{i})*Df{i});
                Pf            {i} = Inv{i}*transpose(Dwf{i});
                pde.mat.cd_Df (i) = cond(Df {i});
                pde.mat.cd_Dwf(i) = cond(Dwf{i});
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
            
            % >> Face quadrature.
            %  > Function handle coordinate transformation.
            [xy_fg,j_fg] = B_1_2.CD_1D();
            for i = 1:msh.f.NF
                pde.f.gq{i} = B_1_2.GaussFace_Points(xy_fg,j_fg,ng,msh.f.xy_v{i});
            end
            
%             %  > Convection.
%             for i = 1:msh.f.NF
%                 df  {i} = SubClass_2_2.Compute_df(1,Face(i,:),msh.f.gq{i},Coeff_1,Exp_1);
%                 Tf_C{i} = df{i}*Pf{i};
%             end
%             %  > Diffusion.
%             for i = 1:msh.f.NF
%                 grad_df{i} = SubClass_2_2.Compute_df(2,Face(i,:),msh.f.gq{i},Coeff_2,Exp_2);
%                 Tf_D   {i} = grad_df{i}*Pf{i};
%             end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        % >> 2.2. ---------------------------------------------------------
    end
end
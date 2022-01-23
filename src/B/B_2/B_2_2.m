classdef B_2_2
    methods (Static)
        %% > Wrap-up B_2_2.
        function [pde] = WrapUp_B_2_2(msh,pde,ng,vx,vy,gx,gy,bnd_ff,bnd_fc)
            % >> ----------------------------------------------------------
            % >> 1.   Assemble matrices.
            %  > 1.1. Assemble matrices Df, Dwf, Pf and Tf.
            %  > 1.2. Compute df array.
            %  > 1.3. Assemble matrices A and B.
            %  > 1.4. Solver setup.
            % >> ----------------------------------------------------------
            % >> 1.
            %  > 1.1.
            pde = B_2_2.Assemble_Mat_1(msh,pde,ng,bnd_ff,bnd_fc);
            %  > 1.3.
            pde = B_2_2.Assemble_Mat_2(msh,pde,vx,vy,gx,gy,bnd_ff);
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [pde] = Assemble_Mat_1(msh,pde,ng,bnd_ff,bnd_fc)
            %% > Matrices Df and Dwf.
            % >> Df = [1*{(x-xf)^0}*{(y-yf)^0},1*{(x-xf)^1}*{(y-yf)^0},1*{(x-xf)^0}*{(y-yf)^1},...] = [1,(x-xf),(y-yf),...].           
            % >> (x,y) coordinates...
            i          = 1:msh.f.NF;
            %  > ...face centroid.
            Face(1,i)  = msh.f.mean(1,i);
            Face(2,i)  = msh.f.mean(2,i);
            %  > ...stencil points.
            stl_xy     = msh.s.xy_v_t;
            
            for i = 1:msh.f.NF
                %  > j.
                j{i} = 1:size(stl_xy{i},2);  
                %  > x-xf.
                XY{i}(1,j{i}) = stl_xy{i}(1,j{i})-Face(1,i);
                XY{i}(2,j{i}) = stl_xy{i}(2,j{i})-Face(2,i);
                %  > Df.
                for k = 1:pde.pr.numb
                    Df{i}(j{i},k) = pde.pr.Coeff_1(k).*(XY{i}(1,j{i}).^pde.pr.Exp_1(1,k)).*(XY{i}(2,j{i}).^pde.pr.Exp_1(2,k));
                end
            end
     
            % >> Dwf = w*Df.
            for i = 1:msh.f.NF
                for j = 1:size(stl_xy{i},2)
                    if ~isequal(stl_xy{i}(:,j),Face(:,i))
                        %  > Point.
                        Point{i}(:,j) = stl_xy{i}(:,j);
                        %  > Distance.
                        d    {i}(j)   = A_Tools.fft_dist_1(Point{i}(:,j),Face(:,i));
                        w    {i}(j)   = pde.wf.wf_1(d{i}(j),1);
                    else
                        %  > Point.
                        %    Remark: If point=face, use cell centroid instead.
                        Point{i}(:,j) = msh.c.mean(:,bnd_fc(bnd_ff == i));
                        %  > Distance.
                        d    {i}(j)   = A_Tools.fft_dist_1(Point{i}(:,j),Face(:,i));
                        w    {i}(j)   = pde.wf.wf_2(d{i}(j),1);
                    end
                    Dwf{i}(j,:) = Df{i}(j,:).*w{i}(j);
                end
            end
            
            %% > Matrices Pf and Tf.
            % >> Pf.
            %  > Pf = inv(Dwf_T*Df)*Dwf_T.
            for i = 1:msh.f.NF
                Inv{i} = pinv(transpose(Dwf{i})*Df{i});
                Pf {i} = Inv{i}*transpose(Dwf{i});
            end
            % >> Tf = [1,(x-xf),(y-yf),...]*Pf*Phi = df*Pf*Phi.
            for i = 1:msh.f.NF
                %  > Face quadrature.
                pde.f.gq{i} = B_1_2.GaussFace_Points(pde.f.xy_fg,pde.f.j_fg,ng,msh.f.xy_v{i});
                %  > Convection.
                df      {i} = B_2_2.Compute_df(1,Face(:,i),pde.f.gq{i},pde.pr.numb,pde.pr.Coeff_1,pde.pr.Exp_1);
                Tf_C    {i} = df{i}*Pf{i};
                %  > Diffusion.
                grad_df {i} = B_2_2.Compute_df(2,Face(:,i),pde.f.gq{i},pde.pr.numb,pde.pr.Coeff_2,pde.pr.Exp_2);
                Tf_D    {i} = grad_df{i}*Pf{i};
            end
            
            % >> Deal...
            %  > Df, Dwf and Pf
            [pde.f.Df,pde.f.Dwf,pde.f.Pf] = ...
                deal(Df,Dwf,Pf);
            %  > Tf_C and .
            [pde.f.Tf_C,pde.f.Tf_D] = ...
                deal(Tf_C,Tf_D);
        end
        % >> 1.2. ---------------------------------------------------------
        function [df_ij] = Compute_df(iD,mean_f,fg,numb,Coeff,Exp)
            if iD == 1
                % >> Phi_f.
                i          = 1:size(fg.Points,1);
                j          = 1:numb;
                df_ij(i,j) = Coeff(j).*fg.Weights(i).*((mean_f(1)-fg.Points(i,1)).^Exp(1,j)).*((mean_f(2)-fg.Points(i,2)).^Exp(2,j));
                df_ij      = sum(df_ij,1);
            elseif iD == 2
                % >> grad(Phi_f).
                for i = 1:size(Coeff,1)
                    j              = 1:size(fg.Points,1);
                    k              = 1:numb;
                    df_ijk{i}(j,k) = Coeff(i,k).*fg.Weights(j).*((mean_f(1)-fg.Points(j,1)).^Exp{i}(1,k)).*((mean_f(2)-fg.Points(j,2)).^Exp{i}(2,k));
                    df_ijk{i}      = sum(df_ijk{i},1);
                    df_ij (i,:)    = [df_ijk{i}];
                end
            end
        end
        % >> 1.3. ---------------------------------------------------------
        function [pde] = Assemble_Mat_2(msh,pde,vx,vy,gx,gy,bnd_ff)            
            % >> X*Tf.
            %    Remark: Tf=[A,B,C,D,E,F,G,H,I,J,...], where: Cell dependent coefficients: A,B,C,...,G.
            %                                                 Face dependent coefficients: H,i,J,...,(...).
            for i = 1:msh.f.NF
                %  > V*Tf_C (face 'i').
                V_Tf_C{i}(1,:) = vx.*pde.f.Tf_C{i};
                V_Tf_C{i}(2,:) = vy.*pde.f.Tf_C{i};
                %  > G*Tf_D (face 'i').
                G_Tf_D{i}(1,:) = gx.*pde.f.Tf_D{i}(1,:);
                G_Tf_D{i}(2,:) = gy.*pde.f.Tf_D{i}(2,:);
                %  > Phi_c.
                Phi_fc{i}      = [msh.s.c{:,i}];
                %  > Phi_f.
                if ~isempty(msh.s.xy_v_f{i})
                    Phi_ff{i}  = [msh.s.f{:,i}];
                end
            end
            % >> Sentil elements.
            %  > X*Tf*Sf.
            for i = 1:msh.c.NC
                for j = 1:size(msh.c.f.Sf{i},2)
                    face     {i}(j) = msh.c.f.f {i}(j);
                    V_Tf_C_Sf{i}{j} = msh.c.f.Sf{i}(:,j)'*V_Tf_C{face{i}(j)};
                    G_Tf_D_Sf{i}{j} = msh.c.f.Sf{i}(:,j)'*G_Tf_D{face{i}(j)};
                end
            end
            
            % >> A (Cell dependent coefficients).
            %  > Initialize.
            A = zeros(msh.c.NC,msh.c.NC);
            for i = 1:msh.c.NC
                %  > ... for cell #i: Phi_f(j_1)=A*Phi(X)+B*Phi(Y)+C*Phi(Z)+...
                %                     Phi_f(j_2)=D*Phi(M)+E*Phi(N)+F*Phi(O)+...
                %                     Phi_f(...)=...
                for j = 1:length(face{i})
                    %  > Number of cells used in the stencil.
                    numb_c{i}(j) = length(Phi_fc{face{i}(j)});
                    %  > (k_ini,k_fin).
                    k_ini {i}(j) = 1;
                    k_fin {i}(j) = numb_c{i}(j);
                    
                    for k = k_ini{i}(j):k_fin{i}(j)
                        %  > Index in the j-direction.
                        ijk = Phi_fc{face{i}(j)}(k);
                        
                        % >> Convection-diffusion.
                        %  > Convection.
                        A(i,ijk) = A(i,ijk)+V_Tf_C_Sf{i}{j}(k);
                        %  > Diffusion.
                        A(i,ijk) = A(i,ijk)-G_Tf_D_Sf{i}{j}(k);
                    end
                end
            end
            A = sparse(A);
            
            % >> B (Face dependent coefficients).
            %  > Initialize.
            B = zeros(msh.c.NC,1);
            for i = 1:msh.c.NC
                for j = 1:length(face{i})
                    %  > Number of faces used in the stencil.
                    numb_f{i}(j) = length(Phi_ff{face{i}(j)});
                    %  > (k_ini,k_fin).
                    k_ini {i}(j) = numb_c{i}(j)+1;
                    k_fin {i}(j) = numb_c{i}(j)+numb_f{i}(j);
                    
                    for k = k_ini{i}(j):k_fin{i}(j)
                        % >> Convection-diffusion.
                        %  > Convection.
                        B(i,1) = B(i,1)+0;
                        %  > Diffusion.
                        B(i,1) = B(i,1)+0;
                    end
                end
            end
            %  > Add source term.
            for i = 1:msh.c.NC
                B(i,1) = B(i,1)+pde.c.F_Vol(i);
            end
            B = sparse(B);
            
            % >> Solve: AX=B.
            pde.Phi = B_2_2.SetUp_bicgstabl(A,B,10e-12,10e3)';
        end
        % >> 1.4. ---------------------------------------------------------
        function [X] = SetUp_bicgstabl(A,B,Tol,iterMax)
            setup = struct('type','ilutp','droptol',1e-6);
            [L,U] = ilu(A,setup);
            [X,~] = bicgstabl(A,B,Tol,iterMax,L,U,[]);
        end
    end
end
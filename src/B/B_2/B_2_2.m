classdef B_2_2
    methods (Static)
        %% > Wrap-up B_2_2.
        function [pde] = WrapUp_B_2_2(msh,pde,vx,vy,gx,gy,bnd_ff,bnd_fc)
            % >> ----------------------------------------------------------
            % >> 1.   Assemble matrices.
            %  > 1.1. Assemble matrices Df, Dwf, Pf and Tf.
            %  > 1.2. Compute df array.
            %  > 1.3. Assemble matrices A and B.
            %  > 1.4. Solver setup.
            % >> ----------------------------------------------------------
            % >> 1.
            %  > 1.1.
            pde = B_2_2.Assemble_Mat_1(msh,pde,bnd_ff,bnd_fc);
            %  > 1.2.
            pde = B_2_2.Assemble_Mat_2(msh,pde,bnd_ff,vx,vy,gx,gy);
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [pde] = Assemble_Mat_1(msh,pde,bnd_ff,bnd_fc)
            %% > Matrices Df and Dwf.
            % >> (x,y) coordinates...
            i         = 1:msh.f.NF;
            %  > ...face centroid.
            Face(1,i) = msh.f.mean(1,i);
            Face(2,i) = msh.f.mean(2,i);
            %  > ...stencil points.
            stl_xy    = msh.s.xy_v_t;
            
            % >> Df.
            for i = 1:msh.f.NF
                j {i} = 1:size(stl_xy{i},2);
                %  > x-xf.
                XY{i}(1,j{i}) = stl_xy{i}(1,j{i})-Face(1,i);
                XY{i}(2,j{i}) = stl_xy{i}(2,j{i})-Face(2,i);
                %  > Df = [1*{(x-xf)^0}*{(y-yf)^0},1*{(x-xf)^1}*{(y-yf)^0},1*{(x-xf)^0}*{(y-yf)^1},...] = [1,(x-xf),(y-yf),...].
                for k = 1:pde.pr.numb
                    Df{i}(j{i},k) = pde.pr.Coeff_1(k).*(XY{i}(1,j{i}).^pde.pr.Exp_1(1,k)).*(XY{i}(2,j{i}).^pde.pr.Exp_1(2,k));
                end
            end
            
            % >> Dwf = w*Df.
            for i = 1:msh.f.NF
                %  > if d~=0...
                for j = 1:size(stl_xy{i},2)
                    d{i}  (j) = A_Tools.fft_dist_1([stl_xy{i}(:,j)';Face(:,i)']);
                end
                w    {i}(:,1) = pde.wf.wf_1(d{i},2);
                %  > if d=0...
                %    Remark: If point=face, use cell centroid instead.
                d_flag{i} = find(~d{i});
                if ~isempty(d_flag{i})
                    k         = d_flag{i};
                    Point     = msh.c.mean(:,bnd_fc(bnd_ff == i));
                    d{i}  (k) = A_Tools.fft_dist_1([Point';Face(:,i)']);
                    w{i}(k,1) = pde.wf.wf_1(d{i}(k),2);
                end
                Dwf{i} = bsxfun(@times,Df{i},w{i});
            end
            
            %% > Matrices Pf and Tf.
            % >> Pf.
            %  > Pf = inv(Dwf_T*Df)*Dwf_T.
            for i = 1:msh.f.NF
                Pf{i} = pinv(transpose(Dwf{i})*Df{i})*transpose(Dwf{i});
            end
            % >> Tf = [1,(x-xf),(y-yf),...]*Pf*Phi = df*Pf*Phi.
            for i = 1:msh.f.NF
                %  > Face quadrature.
                pde.f.gq  {i} = B_1_2.GaussFace_Points(pde.f.Q_1D,pde.f.xy_fg,pde.f.j_fg,msh.f.xy_v{i});
                %  > Convection/diffusion contribution(s).
                df        {i} = B_2_2.Compute_df(1,Face(:,i),pde.f.gq{i},pde.pr.numb,pde.pr.Coeff_1,pde.pr.Exp_1);
                grad_df   {i} = B_2_2.Compute_df(2,Face(:,i),pde.f.gq{i},pde.pr.numb,pde.pr.Coeff_2,pde.pr.Exp_2);
                pde.f.Tf_C{i} = df{i}*Pf{i};
                pde.f.Tf_D{i} = grad_df{i}*Pf{i};
            end
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
        function [pde] = Assemble_Mat_2(msh,pde,bnd_ff,vx,vy,gx,gy)
            %  > Initialize.
            A = zeros(msh.c.NC,msh.c.NC);
            B = zeros(msh.c.NC,1);
            
            % >> X*Tf.
            %    Remark: Tf=[A,B,C,D,E,F,G,H,I,J,...], where: Cell dependent coefficients: A,B,C,...,G.
            %                                                 Face dependent coefficients: H,i,J,...,(...).
            % >> Face 'i'...
            for i = 1:msh.f.NF
                %  > Convection/diffusion contribution(s).
                V_Tf_C{i}(1,:) = vx.*pde.f.Tf_C{i};
                V_Tf_C{i}(2,:) = vy.*pde.f.Tf_C{i};
                G_Tf_D{i}(1,:) = gx.*pde.f.Tf_D{i}(1,:);
                G_Tf_D{i}(2,:) = gy.*pde.f.Tf_D{i}(2,:);
                %  > Cell/Face indices.
                Phi_fc{i}      = [msh.s.c{:,i}];
                if ~isempty(msh.s.xy_v_f{i})
                    Phi_ff{i}  = [msh.s.f{:,i}];
                end
            end
            
            for i = 1:msh.c.NC
                for j = 1:size(msh.c.f.Sf{i},2)
                    %  > Face 'j' of cell 'i'.
                    face     {i}(j) = msh.c.f.f {i}(j);
                    %  > Convection/diffusion contribution(s).
                    V_Tf_C_Sf{i}{j} = msh.c.f.Sf{i}(:,j)'*V_Tf_C{face{i}(j)};
                    G_Tf_D_Sf{i}{j} = msh.c.f.Sf{i}(:,j)'*G_Tf_D{face{i}(j)};
                end
            end
            
            % >> A (Cell dependent coefficients).
            %  > ... for cell #i: Phi_f(j_1)=A*Phi(X)+B*Phi(Y)+C*Phi(Z)+...
            %                     Phi_f(j_2)=D*Phi(M)+E*Phi(N)+F*Phi(O)+...
            %                     Phi_f(...)=...
            for i = 1:msh.c.NC
                for j = 1:length(face{i})
                    %  > Index in the j-direction.
                    numb_c   {i}(j)     = length(Phi_fc{face{i}(j)});
                    kc_i     {i}(j)     = 1;
                    kc_f     {i}(j)     = numb_c{i}(j);
                    k                   = kc_i{i}(j):kc_f{i}(j);
                    ijk_c    {i}{j}     = Phi_fc{face{i}(j)}(k);
                    %  > Convection/diffusion contribution(s).
                    A(i,ijk_c{i}{j}(k)) = A(i,ijk_c{i}{j}(k))+V_Tf_C_Sf{i}{j}(k);
                    A(i,ijk_c{i}{j}(k)) = A(i,ijk_c{i}{j}(k))-G_Tf_D_Sf{i}{j}(k);
                end
            end
            
            % >> B (Face dependent coefficients).
            for i = 1:msh.c.NC
                %  > Boundary contribution(s).
                for j = 1:length(face{i})
                    %  > Index in the j-direction.
                    numb_f   {i}(j) = length(Phi_ff{face{i}(j)});
                    kf_i     {i}(j) = numb_c{i}(j)+1;
                    kf_f     {i}(j) = numb_c{i}(j)+numb_f{i}(j);
                    k               = kf_i{i}(j):kf_f{i}(j);
                    ijk_f    {i}{j} = Phi_ff{face{i}(j)}(k-numb_c{i}(j));
                    %  > Convection/diffusion contribution(s).
                    if all(isempty(ijk_f{i}{j}))
                        %  > Skip cells that don't have faces whose sentil does not include boundary faces.
                        break;
                    else
                        if ~isempty(ijk_f{i}{j})
                            for l = 1:length(ijk_f{i}{j})
                                %  > Indices.
                                bnd_j{i}{j}(l) = find(ijk_f{i}{j}(l) == bnd_ff);
                                %  > Boundary values.
                                bnd_f{i}{j}(l) = pde.bnd.f(bnd_j{i}{j}(l));
                            end
                        end
                    end
                    B(i,1) = B(i,1)-V_Tf_C_Sf{i}{j}(k)*bnd_f{i}{j}';
                    B(i,1) = B(i,1)+G_Tf_D_Sf{i}{j}(k)*bnd_f{i}{j}';
                end
                %  > Source term contribution(s).
                B(i,1) = B(i,1)+pde.c.F_Vol(i);
            end
            
            % >> PDE solution (nodal values).
            pde.Phi = B_2_2.SetUp_bicgstabl(A,B,10e-12,10e3)';
        end
        % >> 1.4. ---------------------------------------------------------
        function [X] = SetUp_bicgstabl(A,B,Tol,iterMax)
            [A,B,setup] = deal(sparse(A),sparse(B),struct('type','ilutp','droptol',1e-6));
            [L,U]       = ilu(A,setup);
            [X,~]       = bicgstabl(A,B,Tol,iterMax,L,U,[]);
        end
    end
end
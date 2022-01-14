classdef SubClass_2_2
    methods (Static)
        %% > Wrap up SubClass_2_2.
        function [msh] = WrapUp_2_2(inp,msh,fn)
            % >> ----------------------------------------------------------
            % >> 1.   Set stencil coordinates (gather face/cell centroids' coordinates).
            %  > 1.1. Wrap up '1.' and call '1.2.'
            %  > 1.2. Find neighbouring faces' centroid coordinates.
            %  > 1.3. Gather stencil elements' coordinates. 
            % >> 2.   Assemble matrices Df, Dwf, Pf and Tf.
            %  > 2.1. Wrap up '2.' and call '2.2.'.
            %  > 2.2. Compute polynomial terms.
            %  > 2.3. Set weighting function.
            %  > 2.4. Quadrature abcissas/weights.
            %  > 2.5. Coordinate transformation: csi->x/y.
            %  > 2.6. Compute df array.
            % >> 3.   Assemble matrices A and B.
            %  > 3.1. Assemble matrix A.
            %  > 3.2. Assemble matrix B.
            % >> ----------------------------------------------------------
            % >> Local variables.
            wf    = inp.fr.wf;
            ng    = inp.fr.ng;
            Ext_1 = inp.fr.Ext_1;
            V     = [inp.pr.vx,inp.pr.vy];
            G     = [inp.pr.gx,inp.pr.gy];
                        
            % >> 1.
            [msh,st_xy] = SubClass_2_2.WrapUp_1(msh,Ext_1);
            % >> 2.
            [msh,Tf_C,Tf_D] = SubClass_2_2.WrapUp_2(msh,st_xy,wf,ng);
            % >> 3.
            [msh] = SubClass_2_2.WrapUp_3(msh,V,G,Tf_C,Tf_D);
        end
        
        %% > Tools.
        %% > 1.) ----------------------------------------------------------
        % >> 1.1.) --------------------------------------------------------
        function [msh,st_xy_i] = WrapUp_1(msh,Ext_1)
            % >> iD's.
            %  > Boundary cells' index.
            for i = 1:size(msh.bnd.c,2)
                bnd_cc(i) = msh.bnd.c{2,i};
            end
            %  > Boundary faces' cell index.
            for i = 1:size(msh.bnd.f,2)
                bnd_ff(i) = msh.bnd.f{2,i};
                bnd_fc(i) = msh.bnd.f{3,i};
            end
            
            % >> Deal (Xv,Yv) coordinates of stencil elements.
            for i = 1:msh.f.NF
                for j = 1:size(msh.s.st_i,1)
                    %  > Stencil cell centroids.
                    len_c{i}(j) = length(msh.s.st_i{j,i});
                    for k = 1:len_c{i}(j)
                        st_v{j,i}(1,k) = msh.c.mean(1,msh.s.st_i{j,i}(k));
                        st_v{j,i}(2,k) = msh.c.mean(2,msh.s.st_i{j,i}(k));
                    end
                    %  > Check whether face i's stencil cells belong to the boundary. If so, add the respective face to the stencil.
                    Flag{j,i} = zeros(1,len_c{i}(j));
                    Flag{j,i} = ismembc(msh.s.st_i{j,i},bnd_cc);
                    if any(Flag{j,i} == 1)
                        %  > Call 1.2.)
                        st_v{j,i} = SubClass_2_2.Deal_FaceCoord(Flag{j,i},st_v{j,i},bnd_ff,bnd_fc,len_c{i}(j),msh.s.st_i{j,i},msh.f.mean);
                    end
                    msh.s.xy_v{j,i} = st_v{j,i};
                end
            end
            %  > Gather stencil elements' coordinates.
            st_xy_i = SubClass_2_2.Gather_Elements(msh);
            
%             % >> > Perform stencil extension?
%             %  > Remark: Use functions on SubClass_1_3.
%             msh = SubClass_1_3.Set_Limits(msh,st_xy_i);
%             
%             %  > Extension #1.
%             if strcmpi(Ext_1,'T')
%                 %  > Elements' index to be added to the stencil.
%                 add_to = SubClass_1_3.StencilExt_1(msh);
%                 %  > Add elements' coordinates...
%                 for i = 1:msh.f.NF
%                     if ~isempty(add_to{i})
%                         msh.s.st_f{i} = add_to{i};
%                         st_xy_i   {i} = [st_xy_i{i},msh.c.mean(:,msh.s.st_f{i})];
%                     end
%                 end
%             end
            
            
            
            
        end
        % >> 1.2.) --------------------------------------------------------
        function [st_v] = Deal_FaceCoord(Flag,st_v,bnd_ff,bnd_fc,len_c,st,mean_f)
            %  > Add respective boundary faces.
            %    Remark: A given boundary cell may contain more than 1 boundary face (see 3rd row of msh.bnd.f)
            j = 1;
            for i = 1:length(Flag)
                if Flag(i)
                    cf{j} = find(bnd_fc == st(i));
                    j     = j+1;
                end
            end
            %  > Faces to be added to the stencil.
            f_add = cell2mat(cf);
            for j = 1:length(f_add)
                st_v(1,j+len_c) = mean_f(1,bnd_ff(f_add(j)));
                st_v(2,j+len_c) = mean_f(2,bnd_ff(f_add(j)));
            end
        end
        % >> 1.3.) --------------------------------------------------------
        function [st_xy] = Gather_Elements(msh)
            for i = 1:size(msh.s.st_i,2)
                for j = 1:size(msh.s.st_i,1)
                    st_xy{i}{j} = msh.s.xy_v{j,i};
                end
                st_xy{i} = cell2mat(st_xy{i});
            end
        end
                
        %% > 2.) ----------------------------------------------------------
        % >> 2.1.) --------------------------------------------------------
        function [msh,Tf_C,Tf_D] = WrapUp_2(msh,st_xy,wf,ng)
            %% > Polynomial reconstruction.
            %  > Polynomial order.
            p = 2.*size(msh.s.st_i,1)-1;
            %  > Number of terms.
            len_p = 1./2.*(p+1).*(p+2);
            %  > Polynomial terms.
            [Coeff_1,Exp_1] = SubClass_2_2.Compute_PolTerms(1,p,len_p);
            [Coeff_2,Exp_2] = SubClass_2_2.Compute_PolTerms(2,p,len_p);
                 
            %% > Matrices Df and Dwf.
            % >> Df.
            %  > Df = [1*{(x-xf)^0}*{(y-yf)^0},1*{(x-xf)^1}*{(y-yf)^0},1*{(x-xf)^0}*{(y-yf)^1},...] = [1,(x-xf),(y-yf),...].
            for i = 1:msh.f.NF
                Face(i,:) = ...
                    [msh.f.mean(1,i),msh.f.mean(2,i)];
                for j = 1:size(st_xy{i},2)
                    XY{i}(1,j)  = st_xy{i}(1,j)-Face(i,1);
                    XY{i}(2,j)  = st_xy{i}(2,j)-Face(i,2);
                    for k = 1:len_p
                        Df{i}(j,k) = Coeff_1(k).*(XY{i}(1,j).^Exp_1(1,k)).*(XY{i}(2,j).^Exp_1(2,k));
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
                for j = 1:size(st_xy{i},2)
                    %  > Point.
                    %  > Remark: If point=face, then face "i" belongs to a corner cell -> use cell centroid instead.
                    if ~isequal(st_xy{i}(:,j)',Face(i,:))
                        Point{i}(j,:) = st_xy{i}(:,j)';
                        w    {i}(j,1) = SubClass_2_2.W_Function(Point{i}(j,:),Face(i,:),wf,4,'F');
                    else
                        Point{i}(j,:) = msh.c.mean(:,bnd_fc(bnd_ff == i));
                        w    {i}(j,1) = SubClass_2_2.W_Function(Point{i}(j,:),Face(i,:),wf,4,'T');
                    end
                    Dwf{i}(j,:) = Df{i}(j,:).*w{i}(j,1);
                end
            end
            
            %% > Matrices Pf and Tf. 
            % >> Pf.
            %  > Pf = inv(Dwf_T*Df)*Dwf_T.
            for i = 1:msh.f.NF
                Inv          {i} = eMatrices(transpose(Dwf{i})*Df{i});
                Pf           {i} = Inv{i}*transpose(Dwf{i});
                msh.s.cond_Df(i) = cond(Df{i});
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
            [xy_fg,j_fg] = SubClass_2_2.Create_FunctionHandle();
            for i = 1:msh.f.NF
                msh.f.gq{i} = SubClass_2_2.GaussFace_Points(xy_fg,j_fg,ng,msh.f.xy_v{i});
            end
            
            %  > Convection.
            for i = 1:msh.f.NF
                df  {i} = SubClass_2_2.Compute_df(1,Face(i,:),msh.f.gq{i},Coeff_1,Exp_1);
                Tf_C{i} = df{i}*Pf{i};
            end
            %  > Diffusion.
            for i = 1:msh.f.NF
                grad_df{i} = SubClass_2_2.Compute_df(2,Face(i,:),msh.f.gq{i},Coeff_2,Exp_2);
                Tf_D   {i} = grad_df{i}*Pf{i};
            end
        end
        % >> 2.2.) --------------------------------------------------------
        function [Coeff_iD,Exp_iD] = Compute_PolTerms(iD,p,len_p)
            %                |------------------------------------------------------------------------------------------------|
            %  > Term      = | 1 | x | y | x^2 | xy  | y^2 | x^3  | x^2y | xy^2 | y^3  | x^4  | x^3y  | x^2y^2 | xy^3  | y^4  | (...)
            %                |------------------------------------------------------------------------------------------------|
            %  > Coeff     = | 1 | 1 | 1 | 1   | 1   | 1   | 1    | 1    | 1    | 1    | 1    | 1     | 1      | 1     | 1    | (...)
            %  > X_exp     = | 0 | 1 | 0 | 2   | 1   | 0   | 3    | 2    | 1    | 0    | 4    | 3     | 2      | 1     | 0    | (...)
            %  > Y_exp     = | 0 | 0 | 1 | 0   | 1   | 2   | 0    | 1    | 2    | 3    | 0    | 1     | 2      | 3     | 4    | (...)
            %                |------------------------------------------------------------------------------------------------|
            %                |------------------------------------------------------------------------------------------------|
            %  > dX        = | 0 | 1 | 0 | 2x  | y   | 0   | 3x^2 | 2xy  | y^2  | 0    | 4x^3 | 3x^2y | 2xy^2  | y^3   | 0    | (...)
            %  > Coeff     = | 0 | 1 | 0 | 2   | 1   | 0   | 3    | 2    | 1    | 0    | 4    | 3     | 2      | 1     | 0    | (...)
            %  > dX_exp(x) = | 0 | 0 | 0 | 1   | 0   | 0   | 2    | 1    | 0    | 0    | 3    | 2     | 1      | 0     | 0    | (...)
            %  > dX_exp(y) = | 0 | 0 | 0 | 0   | 1   | 0   | 0    | 1    | 2    | 0    | 0    | 1     | 2      | 3     | 0    | (...)            
            %                |------------------------------------------------------------------------------------------------|
            %  > dY        = | 0 | 0 | 1 | 0   | x   | 2y  | 0    | x^2  | 2xy  | 3y^2 | 0    | x^3   | 2x^2y  | 3xy^2 | 4y^3 | (...)
            %  > Coeff     = | 0 | 0 | 1 | 0   | 1   | 2   | 0    | 1    | 2    | 3    | 0    | 1     | 2      | 3     | 4    | (...)
            %  > dX_exp(x) = | 0 | 0 | 0 | 0   | 1   | 0   | 0    | 2    | 1    | 0    | 0    | 3     | 2      | 1     | 0    | (...)
            %  > dX_exp(y) = | 0 | 0 | 0 | 0   | 0   | 1   | 0    | 0    | 1    | 2    | 0    | 0     | 1      | 2     | 3    | (...)
            %                |------------------------------------------------------------------------------------------------|
            %                |------------------------------------------------------------------------------------------------|      
            %  > #         = | 1 | 2 | 3 | 4   | 5   | 6   | 7    | 8    | 9    | 10   | 11   | 12    | 13     | 14    | 15   | (...)
            %                |------------------------------------------------------------------------------------------------| (...)
            %                | p = 1     | p = 3                                | p = 5                                      ...
            %                |------------------------------------------------------------------------------------------------|
            
            % >> iD = 1.
            %  > Polynomial coefficients.
            Coeff_1 = ones(1,len_p);
            %  > Polynomial exponents.            
            Exp_1   = zeros(2,len_p);
            i       = 1;
            while i < p+1
                %  > Previous last index.
                k = 1./2.*(i.^2+i);
                %  > X.
                for j = 1:i+1
                    Exp_1(1,k+j) = i+1-j;
                end
                %  > Y.
                for j = 1:i+1
                    Exp_1(2,k+j) = j-1;
                end
                i = i+1;
            end 
            
            % >> iD = 2.
            %  > Polynomial exponents.
            if iD == 2
                [Exp_dX,Exp_dY] = deal(zeros(2,len_p));
                %  > grad_x(x): Polynomial exponents.
                i = 2;
                while i < p+1
                    %  > Previous last index.
                    k = 1./2.*(i.^2+i);
                    %  > X.
                    for j = 1:i
                        Exp_dX(1,k+j) = i-j;
                    end
                    %  > Y.
                    for j = 1:i
                        Exp_dX(2,k+j) = j-1;
                    end
                    i = i+1;
                end
                %  > grad_x(y): Polynomial exponents.
                i = 2;
                while i < p+1
                    %  > Previous last index.
                    k = 1./2.*(i.^2+i);
                    %  > X.
                    for j = 1:i
                        Exp_dY(1,k+j+1) = i-j;
                    end
                    %  > Y.
                    for j = 1:i
                        Exp_dY(2,k+j+1) = j-1;
                    end
                    i = i+1;
                end
            end

            % >> Deal coefficients/exponents.
            if iD == 1
                % >> Phi_f.
                Coeff_iD = Coeff_1;
                Exp_iD   = Exp_1;
            elseif iD == 2
                % >> grad(Phi_f).
                %  > grad(x).
                Coeff_iD (1,:) = Exp_1 (1,:);
                Exp_iD{1}(1,:) = Exp_dX(1,:);
                Exp_iD{1}(2,:) = Exp_dX(2,:);
                %  > grad(y).
                Coeff_iD (2,:) = Exp_1 (2,:);
                Exp_iD{2}(1,:) = Exp_dY(1,:);
                Exp_iD{2}(2,:) = Exp_dY(2,:);
            end
        end       
        % >> 2.3.) --------------------------------------------------------
        function [W_func] = W_Function(Point,Face,wf,k,Flag)
            % >> Weighting functions: 1) Unweighted: wij = 1.
            %                         2) Weighted  : wij = 1/|xj-xi|^2 -> wij = 1/|dist(x_point-xf)|^k.
            %                         3) Weighted  : wij = (...).
            if strcmpi(wf,'Unweighted')
                %  > 1)
                W_func = 1;
            elseif strcmpi(wf,'Weighted')
                %  > 2)
                d = pdist([Point;Face],'euclidean');
                if strcmpi(Flag,'F')
                    W_func = 1./d.^k;
                elseif strcmpi(Flag,'T')
                    W_func = 1./(d./2).^k;
                    A = 1-exp(-k.^2);
                    W_func = 1;
                end
            end
        end
        
        % >> 2.6.) --------------------------------------------------------
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
                       
        %% > 3.) ----------------------------------------------------------
        % >> 3.1.) --------------------------------------------------------
        function [msh] = WrapUp_3(msh,V,G,Tf_C,Tf_D)
            % >> Assemble matrices A and B.
            %  > A.
            SubClass_2_2.Assemble_A(V,G,Tf_C,Tf_D);
            %  > B.
            SubClass_2_2.Assemble_B(V,G,Tf_C,Tf_D);

        end
        % >> 3.2.) --------------------------------------------------------
        function [] = Assemble_A(V,G,Tf_C,Tf_D)
    
        end
        % >> 3.3.) --------------------------------------------------------
        function [] = Assemble_B(V,G,Tf_C,Tf_D)
        end
    end
end
classdef SubClass_2_2
    methods (Static)
        %% > Wrap up SubClass_2_2.
        function [X,Norm] = WrapUp_2_2(inp,msh,fn,bnd,blk)
            %  > Compute source term contribution.
           
            ft = inp.fr.ft;
            st = inp.fr.st;
            np = inp.fr.np;
            ng = inp.fr.ng;
            v = inp.pr.v;
            g = inp.pr.g;
            
            [F_Vol] = SubClass_2_1.Compute_SourceTerm(np,fn,msh);
            
            %  > Set flux reconstruction method.
            if strcmpi(ft,'Explicit')
                [X,Norm] = SubClass_2_2.Reconstruct_ExplicitFlux(msh,bnd,blk,F_Vol);
            elseif strcmpi(ft,'Implicit')
                [X,Norm] = SubClass_2_2.Reconstruct_ImplicitFlux(msh,bnd,blk,F_Vol,np,v,g);
            end
        end
        
        %% > Step #1: Reconstruct faces(Phi,gradPhi).
        % >> NOTE:
        %  > Rather than having multiple files (one for each order), the following procedure is recursively executed based on the method's order (n).
        % >> TODO:
        %  > Solve system: ax=b, where a stores the centroid's distances to the face to be reconstructed and b=[1,0,...]->Phi,b=[0,1,...]->gradPhi.
        %  > c stores the cells used on this reconstruction.
        %
        %  ---------------c----------------
        %    1   2   3        N-2 N-1  N
        %  |-o-|-o-|-o-|-...-|-o-|-o-|-o-|
        %  1   2   3   4    N-2 N-1  N  N+1
        %  ---------------v----------------
        function [a_i,c_i,x] = Reconstruct_Faces(msh,n,iD)
            %% > a(iX).
            % >> Row(s): 1.
            for i = 1:1
                for j = 1:msh.c.NC+1
                    a_i{j}(i,1:n) = 1;
                end
            end
            % >> Row(s): 2,...n.
            for k = 2:n
                %  > Left boundary (downwind).
                for i = 1:n./2
                    j_ini_lhs = -i;
                    j_fin_lhs = -i+n-1;
                    %  > Columns: 2,...,n.
                    for j = j_ini_lhs:j_fin_lhs
                        if j == -i
                            %  > NOTE: Use west boundary "face node".
                            a_i{i}(k,1)     = (msh.f.Xv(1)-msh.f.Xv(i)).^(k-1);
                            c_i{i}(1,1)     = 0;
                        else
                            a_i{i}(k,j+i+1) = (msh.c.Xc(i+j)-msh.f.Xv(i)).^(k-1);
                            c_i{i}(1,j+i+1) = i+j;
                        end
                    end
                end
                %  > Bulk of the domain (central).
                for i = n./2+1:msh.c.NC+1-n./2
                    j_ini_blk = -n./2;
                    j_fin_blk =  n/2-1;
                    %  > Columns: 2,...,n.
                    for j = j_ini_blk:j_fin_blk
                        a_i{i}(k,j+n./2+1) = (msh.c.Xc(i+j)-msh.f.Xv(i)).^(k-1);
                        c_i{i}(1,j+n./2+1) = i+j;
                    end
                end
                %  > Right boundary (upwind).
                for i = msh.c.NC+1-n./2+1:msh.c.NC+1
                    j_ini_rhs = -n+(msh.c.NC+2-i);
                    j_fin_rhs =  msh.c.NC-i+1;
                    %  > Columns: 2,...,n.
                    for j = j_ini_rhs:j_fin_rhs
                        if j == msh.c.NC-i+1
                            %  > NOTE: Use west boundary "face node".
                            a_i{i}(k,n)                = (msh.f.Xv(msh.c.NC+1)-msh.f.Xv(i)).^(k-1);
                            c_i{i}(1,n)                = msh.c.NC+1;  
                        else
                            a_i{i}(k,n-(msh.c.NC+1-i)+j) = (msh.c.Xc(i+j)-msh.f.Xv(i)).^(k-1);
                            c_i{i}(1,n-(msh.c.NC+1-i)+j) = i+j;
                        end
                    end
                end
            end
            %% > b(iX).
            for i = 1:msh.c.NC+1
                b_i{i}       = zeros(n,1);
                b_i{i}(iD,1) = 1;
            end
            %% > x.
            %  > Compute reconstruction coefficients (X).
            %  > 2nd order: x = | B  , A   |           -> Face 1.
            %                   | B  , A   |           -> Face 2.                                                 
            %                   | ..., ... |           -> Face NC+1.
            %  > 4th order: x = | D  , C  , B  , A   | -> Face 1.
            %                   | D  , C  , B  , A   | -> Face 2.                                                 
            %                   | ..., ..., ..., ... | -> Face NC+1.
            %
            %  > NOTE: In order to properly assemble matrix 'a', it is important to be aware of the nodes used for reconstructing the faces.
            %
            %                   | 0  , 1   |           -> Nodes used for reconstruction.
            %  > 2nd order:     | B  , A   |           -> Face 1.
            %                   | 1  , 2   |           -> Nodes used for reconstruction.
            %                   | B  , A   |           -> Face 2.
            %
            %                   | 0  , 1  , 2  , 3   | -> Nodes used for reconstruction.
            %  > 4th order:     | D  , C  , B  , A   | -> Face 1 and 2.
            %                   | 1  , 2  , 3  , 4   | -> Nodes used for reconstruction.
            %                   | D  , C  , B  , A   | -> Face 2.
            for i = 1:msh.c.NC+1
                y{i}   = [a_i{i},b_i{i}];
                x(i,:) = Gauss_Jordan(a_i{i},b_i{i});
            end
        end
        
        %% > Step #2: Assemble matrix a.
        % >> NOTE: 
        %  > Boundary nodes are reconstructed with Phi(0) and Phi(NC+1)! Be aware of the coefficients used when assembling matrix 'a'.
        %  > For example, for cell iX=1 of a 4th order scheme, the LHS and RHS faces are reconstructed with the same nodes. However, for cell iX=2, the RHS face is reconstructed with the RHS nodes of the LHS face.  
        %
        % >> TODO:
        %  > Assemble matrix a.
        %
        % >> for bulk nodes...
        % >> 2nd order: dPhi(i)=Phi(Cell_i,Face_e)-Phi(Cell_i,Face_w)=[A(i+1)*Phi(i+1)+B(i+1)*Phi(i)]-[A(i)*Phi(i)+B(i)*Phi(i-1)] =
        %  > 2nd order:                                              =[Phi(i-1)] *(-B(i)) +[Phi(i)] *(B(i+1) -A(i))+[Phi(i+1)] *(A(i+1))
        %
        %  > a:
        %  > a(i,:)=[...,-B(i),B(i+1)-A(i),A(i+1),...]=[...,i-1,i,i+1,...]              
        % -----------------------------------------------------------------
        % >> 4th order: dPhi(i)=Phi(Cell_i,Face_e)-Phi(Cell_i,Face_w)=[A(i+1)*Phi(i+2)+B(i+1)*Phi(i+1)+C(i+1)*Phi(i)+D(i+1)*Phi(i-1)]-[A(i)*Phi(i+1)+B(i)*Phi(i)+C(i)*Phi(i-1)+D(i)*Phi(i-2)] =
        %                                                            =[Phi(i-2)]*(-D(i))     +[Phi(i-1)]*(D(i+1)-C(i))+[Phi(i)]*(C(i+1)-B(i))+[Phi(i+1)]*(B(i+1)-A(i))+[Phi(i+2)]*(A(i+1))                                  
        %
        %  > a:
        %  > a(i,:)=[...,-D(i),D(i+1)-C(i),C(i+1)-B(i),B(i+1)-A(i),A(i+1),...]=[...,i-2,i-1,i,i+1,i+2,...]
        % -----------------------------------------------------------------
        % >> 6th order: dPhi(i)=Phi(Cell_i,Face_e)-Phi(Cell_i,Face_w)=[A(i+1)*Phi(i+3)+B(i+1)*Phi(i+2)+C(i+1)*Phi(i+1)+D(i+1)*Phi(i)+E(i+1)*Phi(i-1)+F(i+1)*Phi(i-2)]-[A(i)*Phi(i+2)+B(i)*Phi(i+1)+C(i)*Phi(i)+D(i)*Phi(i-1)+E(i)*Phi(i-2)+F(i)*Phi(i-3)] =
        %                                                            =[Phi(i-3)]*(-F(i))+[Phi(i-2)]*(F(i+1)-E(i))+[Phi(i-1)]*(E(i+1)-D(i))+[Phi(i)]*(D(i+1)-C(i))+[Phi(i+1)]*(C(i+1)-B(i))+[Phi(i+2)]*(B(i+1)-A(i))+[Phi(i+3)]*(A(i+1))
        %  > a:
        %  >
        %  > a(i,:)=[...,-F(i),F(i+1)-E(i),E(i+1)-D(i),D(i+1)-C(i),C(i+1)-B(i),B(i+1)-A(i),A(i+1),...]=[...,i-3,i-2,i-1,i,i+1,i+2,i+3,...]
        % -----------------------------------------------------------------
        % >> (...)
        function [a] = Assemble_a(msh,n,xf)
            %% > w/ Bulk nodes ONLY.
            %  > Loop through cells.
            for i = 1:msh.c.NC
                %% > a.
                % >> Fill upper diagonals(i+1,...).
                if i <= n./2-1 || i >= msh.c.NC+1-(n./2-1)
                    %  > w/ boundary nodes.
                    %  > Proceed...
                else
                    for j = 1:n./2
                        if i+j < msh.c.NC+1
                            %  > w/ bulk nodes.
                            if j == n./2
                                a(i,i+j) = xf.Eq(i+1,n);
                            else
                                a(i,i+j) = xf.Eq(i+1,n./2+j)-xf.Eq(i,n./2+j+1);
                            end
                        else
                            %  > Proceed...
                        end
                    end
                end
                % >> Fill main diagonal(i).
                if i <= n./2-1
                    %  > w/ boundary nodes.
                    %  > Proceed...
                elseif i >= msh.c.NC+1-(n./2-1)
                    %  > w/ boundary nodes.
                    %  > Proceed...
                else
                    %  > w/ bulk nodes.
                    a(i,i) = xf.Eq(i+1,n./2)-xf.Eq(i,n./2+1);
                end
                % >> Fill lower diagonals(...,i-1).
                if i <= n./2-1 || i >= msh.c.NC+1-(n./2-1)
                    %  > w/ boundary nodes.
                    %  > Proceed...
                else
                    for j = 1:n./2
                        if i-j > 0
                            %  > w/ bulk nodes.
                            if j == n./2
                                a(i,i-j) = -xf.Eq(i,1);
                            else
                                a(i,i-j) =  xf.Eq(i+1,n./2-j)-xf.Eq(i,n./2-j+1);
                            end
                        else
                            %  > Proceed...
                        end
                    end
                end
            end
            %% > w/ boundary nodes ONLY.
            for i = 1:msh.c.NC
                %% > a.
                % >> Fill upper diagonals(i+1,...).
                if i <= n./2-1 || (i >= msh.c.NC+1-(n./2-1) && i < msh.c.NC)
                    if i <= n./2-1
                        for j = 1:n-(i+1)
                            %  > w/ boundary nodes(EB:i+1).
                            a(i,i+j) = xf.Eq(i+1,i+1+j)-xf.Eq(i,i+1+j);
                        end
                    else
                        for j = msh.c.NC-i:-1:1
                            %  > w/ boundary nodes(WB:i+1).
                            a(i,msh.c.NC+1-j) = xf.Eq(i+1,n-j)-xf.Eq(i,n-j);
                        end
                    end
                else
                    %  > w/ bulk nodes.
                    %  > Proceed...
                end
                % >> Fill main diagonal(i).
                if i <= n./2-1
                    %  > w/ boundary nodes.
                    a(i,i) = xf.Eq(i+1,i+1)-xf.Eq(i,i+1);
                elseif i >= msh.c.NC+1-(n./2-1)
                    %  > w/ boundary nodes.
                    a(i,i) = xf.Eq(i+1,n-1-(msh.c.NC-i))-xf.Eq(i,n-1-(msh.c.NC-i));
                else
                    %  > w/ bulk nodes.
                    %  > Proceed...
                end
                % >> Fill lower diagonals(...,i-1).
                if (i <= n./2-1 && i > 1) || i >= msh.c.NC+1-(n./2-1)
                    if i >= msh.c.NC+1-(n./2-1)
                        for j = 1:n-2-(msh.c.NC-i)
                            %  > w/ boundary nodes.
                            a(i,i-j) = xf.Eq(i+1,n-1-(msh.c.NC-i)-j)-xf.Eq(i,n-1-(msh.c.NC-i)-j);
                        end
                    else
                        for j = 1:i-1
                            %  > w/ boundary nodes.
                            a(i,i-j) = xf.Eq(i+1,i+1-j)-xf.Eq(i,i+1-j);
                        end
                    end
                else
                    %  > w/ bulk nodes.
                    %  > Proceed...
                end
            end
%             % >> Add boundary conditions contribution.
%             %  > West boundary.
%             if strcmpi(obj.West_Boundary,'Neumann')
%                 a = SubClass_2_2.Neumann_BC_a('i_Min',n,a,msh,xf);
%             elseif strcmpi(obj.West_Boundary,'Robin')
%                 a = SubClass_2_2.Robin_BC_a('i_Min',n,a,msh,xf);
%             end
%             %  > East boundary.
%             if strcmpi(obj.East_Boundary,'Neumann')
%                 a = SubClass_2_2.Neumann_BC_a('i_Max',n,a,msh,xf);
%             elseif strcmpi(obj.East_Boundary,'Robin')
%                 a = SubClass_2_2.Robin_BC_a('i_Max',n,a,msh,xf);
%             end
            %  > Transform to sparse notation...
            a = sparse(a);
        end
        
        %% > Step #3: Assemble vector b.
        % >> TODO:
        %  > Assemble vector b.
        function [b] = Assemble_b(msh,n,xf,bnd)
            % >> Initialize b.
            b = zeros(msh.c.NC,1);

            % >> Add boundary conditions contribution.
            %  > West boundary.
%            if strcmpi(obj.West_Boundary,'Dirichlet')
                b = SubClass_2_2.Dirichlet_BC_b('i_Min',n,b,msh,xf,bnd(1));
          %  elseif strcmpi(obj.West_Boundary,'Neumann')
           %     b = SubClass_2_2.Neumann_BC_b('i_Min',n,b,msh,xf,bnd.df(1));
           % elseif strcmpi(obj.West_Boundary,'Robin')
           %     b = SubClass_2_2.Robin_BC_b('i_Min',n,b,msh,xf,bnd.f(1),bnd.df(1));
          %  end
            %  > East boundary.
          %  if strcmpi(obj.East_Boundary,'Dirichlet')
                b = SubClass_2_2.Dirichlet_BC_b('i_Max',n,b,msh,xf,bnd(msh.f.NF));
          %  elseif strcmpi(obj.East_Boundary,'Neumann')
          %      b = SubClass_2_2.Neumann_BC_b('i_Max',n,b,msh,xf,bnd.df(msh.NV));
          %  elseif strcmpi(obj.East_Boundary,'Robin')
          %      b = SubClass_2_2.Robin_BC_b('i_Max',n,b,msh,xf,bnd.f(msh.NV),bnd.df(msh.NV));
          %  end
            b = sparse(b);
        end
        
        %% > Step #4  : Set boundary conditions.     
        %% > Step #4.1: Set a Dirichlet boundary condition (update vector b ONLY).
        % >> Step #4.1: Update vector b.
        function [b] = Dirichlet_BC_b(x_Loc,n,b,msh,xf,Value)
            if strcmpi(x_Loc,'i_Min')
                % >> West boundary ('i_Min').
                for j = 1:n./2
                    if j ~= n./2 && n~=2
                        b(j) = b(j)-(xf.Eq(j+1,1)-xf.Eq(j,1)).*Value;
                    else
                        %  > Last element(LHS).
                        b(j) = b(j)-(-xf.Eq(j,1)).*Value;
                    end
                end
            elseif strcmpi(x_Loc,'i_Max')
                % >> East boundary ('i_Max').
                for j = 1:n./2 
                    if j == 1 || n==2
                        %  > First element(RHS).
                        b(msh.c.NC-n./2+j) = b(msh.c.NC-n./2+j)-(xf.Eq(msh.c.NC-n./2+j+1,n)).*Value;
                    else
                        b(msh.c.NC-n./2+j) = b(msh.c.NC-n./2+j)-(xf.Eq(msh.c.NC-n./2+j+1,n)-xf.Eq(msh.c.NC-n./2+j,n)).*Value;
                    end
                end
            end
        end
        
        %% > Step #4.2: Set a Neumann boundary condition (update matrix a and vector b).
        % >> Step #4.2: Update matrix a.
        function [a] = Neumann_BC_a(x_Loc,n,a,msh,xf)
            if strcmpi(x_Loc,'i_Min')
                % >> West boundary ('i_Min').
                for i = 1:n./2
                    if i ~= n./2 && n ~= 2
                        for j = 1:n-1
                            a(i,j) = a(i,j)+(-xf.gradPhi(1,j+1)./xf.gradPhi(1,1)).*(xf.Eq(i+1,1)-xf.Eq(i,1));
                        end
                    else
                        for j = 1:n-1
                            a(i,j) = a(i,j)+(-xf.gradPhi(1,j+1)./xf.gradPhi(1,1)).*(-xf.Eq(i,1));
                        end
                    end
                end
            elseif strcmpi(x_Loc,'i_Max')
                % >> East boundary ('i_Max').
                for i = msh.c.NC-n./2+1:msh.c.NC
                    if i == msh.c.NC-n./2+1 || n == 2
                        for j = 1:n-1
                            a(i,msh.c.NC-n+1+j) = a(i,msh.c.NC-n+1+j)+(-xf.gradPhi(msh.c.NC+1,j)./xf.gradPhi(msh.c.NC+1,n)).*(xf.Eq(i+1,n));
                        end
                    else
                        for j = 1:n-1
                            a(i,msh.c.NC+1-n+j) = a(i,msh.c.NC+1-n+j)+(-xf.gradPhi(msh.c.NC+1,j)./xf.gradPhi(msh.c.NC+1,n)).*(xf.Eq(i+1,n)-xf.Eq(i,n));
                        end
                    end
                end
            end
        end
        % >> Step #4.2: Update vector b.
        function [b] = Neumann_BC_b(x_Loc,n,b,msh,xf,Value)
            if strcmpi(x_Loc,'i_Min')
                % >> West boundary ('i_Min').
                for i = 1:n./2
                    if i ~= n./2 && n ~= 2
                        b(i) = b(i)-Value./xf.gradPhi(1,1).*(xf.Eq(i+1,1)-xf.Eq(i,1));
                    else
                        b(i) = b(i)-Value./xf.gradPhi(1,1).*(-xf.Eq(i,1));
                    end                    
                end
            elseif strcmpi(x_Loc,'i_Max')
                % >> East boundary ('i_Max').
                for i = msh.c.NC-n./2+1:msh.c.NC
                    if i == msh.c.NC-n./2+1
                        b(i) = b(i)-Value./xf.gradPhi(msh.c.NC+1,n).*(xf.Eq(i+1,n));
                    else
                        b(i) = b(i)-Value./xf.gradPhi(msh.c.NC+1,n).*(xf.Eq(i+1,n)-xf.Eq(i,n));
                    end
                end
            end
        end

        %% > Step #4.3: Set a Robin boundary condition (update matrix a and vector b).
        % >> Step #4.3: Update matrix a.
        function [a] = Robin_BC_a(x_Loc,n,a,msh,xf)
            %  > Deal local variables.
            [V,G] = deal(obj.V,obj.Gamma); 

            if strcmpi(x_Loc,'i_Min')
                % >> West boundary ('i_Min').
                for i = 1:n./2
                    if i ~= n./2 && n~=2
                        for j = 1:n-1
                            a(i,j) = a(i,j)+(G.*xf.gradPhi(1,j+1)./(V-G.*xf.gradPhi(1,1))).*(xf.Eq(i+1,1)-xf.Eq(i,1));
                        end
                    else
                        for j = 1:n-1
                            a(i,j) = a(i,j)+(G.*xf.gradPhi(1,j+1)./(V-G.*xf.gradPhi(1,1))).*(-xf.Eq(i,1));
                        end
                    end
                end
            elseif strcmpi(x_Loc,'i_Max')
                % >> East boundary ('i_Max').
                for i = msh.c.NC-n./2+1:msh.c.NC
                    if i == msh.c.NC-n./2+1 || n == 2
                        for j = 1:n-1
                            a(i,msh.c.NC-n+1+j) = a(i,msh.c.NC-n+1+j)+(G.*xf.gradPhi(msh.c.NC+1,j)./(V-G.*xf.gradPhi(msh.c.NC+1,n))).*(xf.Eq(i+1,n));
                        end
                    else
                        for j = 1:n-1
                            a(i,msh.c.NC+1-n+j) = a(i,msh.c.NC+1-n+j)+(G.*xf.gradPhi(msh.c.NC+1,j)./(V-G.*xf.gradPhi(msh.c.NC+1,n))).*(xf.Eq(i+1,n)-xf.Eq(i,n));
                        end
                    end
                end
            end
        end
        % >> Step #4.3: Update vector b.
        function [b] = Robin_BC_b(x_Loc,n,b,msh,xf,Value_f,Value_df)
            % >> Deal local variables.
            [V,G] = deal(obj.V,obj.Gamma);
            
            if strcmpi(x_Loc,'i_Min')
                % >> West boundary ('i_Min').
                for i = 1:n./2
                    if i ~= n./2 && n~=2
                        b(i) = b(i)-((V.*Value_f-G.*Value_df)./(V-G.*xf.gradPhi(1,1))).*(xf.Eq(i+1,1)-xf.Eq(i,1));
                    else
                        b(i) = b(i)-((V.*Value_f-G.*Value_df)./(V-G.*xf.gradPhi(1,1))).*(-xf.Eq(i,1));
                    end                    
                end
            elseif strcmpi(x_Loc,'i_Max')
                % >> East boundary ('i_Max').
                for i = msh.c.NC-n./2+1:msh.c.NC
                    if i == msh.c.NC-n./2+1
                        b(i) = b(i)-((V.*Value_f-G.*Value_df)./(V-G.*xf.gradPhi(msh.c.NC+1,n))).*(xf.Eq(i+1,n));
                    else
                        b(i) = b(i)-((V.*Value_f-G.*Value_df)./(V-G.*xf.gradPhi(msh.c.NC+1,n))).*(xf.Eq(i+1,n)-xf.Eq(i,n));
                    end
                end
            end
        end
        
        %% > Step #5  : Flux reconstruction.
        %% > Step #5.0: Solver setup.
        function [X] = SetUp_bicgstabl(a,b,Tol,iterMax)
            setup = struct('type','ilutp','droptol',1e-6);
            [L,U] = ilu(a,setup);
            [X,~] = bicgstabl(a,b,Tol,iterMax,L,U,[]);
        end
        %% > Step #5.1: Explicit flux reconstruction.
        % >> Explicit flux reconstruction: Ax-b=0^N.
        function [X,Norm] = Reconstruct_ExplicitFlux(msh,bnd,blk,F_Vol)
            %  > #1: Compute reconstruction coefficients(Phi,gradPhi).
            [~,~,xf.Phi]     = SubClass_2_2.Reconstruct_Faces(msh,obj.n,1);
            [~,~,xf.gradPhi] = SubClass_2_2.Reconstruct_Faces(msh,obj.n,2);
            xf.Eq            = obj.V.*xf.Phi-obj.Gamma.*xf.gradPhi;
            %  > #2: Assemble matrix(a).
            a                = SubClass_2_2.Assemble_a(msh,obj.n,xf);
            b                = SubClass_2_2.Assemble_b(msh,obj.n,xf,bnd);
            %  > #3: Compute approximate solution.
            X.Phi            = blk.f;
            X.Ap             = a*(X.Phi)'-b-F_Vol.Ap';
            X.Ap             = X.Ap';
            %  > #4: Compute error norms.
            X.Error          = abs(zeros(1,msh.c.NC)-X.Ap);
            Norm             = SubClass_2_2.Compute_ErrorNorms(X.Error,msh);
        end
        
        %% > Step #5.2: Implicit flux reconstruction.
        % >> Implicit flux reconstruction: Ax=b, where x=?
        function [X,Norm] = Reconstruct_ImplicitFlux(msh,bnd,blk,F_Vol,np,v,g)
            %  > #1: Compute reconstruction coefficients(Phi,gradPhi).
            [~,~,xf.Phi]     = SubClass_2_2.Reconstruct_Faces(msh,np,1);
            [~,~,xf.gradPhi] = SubClass_2_2.Reconstruct_Faces(msh,np,2);
            xf.Eq            = v.*xf.Phi-g.*xf.gradPhi;
            %  > #2: Assemble matrices(a,b).
            a                = SubClass_2_2.Assemble_a(msh,np,xf);
            b                = SubClass_2_2.Assemble_b(msh,np,xf,bnd);
            b                = b+F_Vol.Ap';
            %  > #3: Compute approximate solution.
            X.Phi            = blk;
            X.Phi_PDE        = SubClass_2_2.SetUp_bicgstabl(a,b,10e-12,10e3);
            %  > #4: Compute error norms.
            X.Error          = abs(X.Phi-X.Phi_PDE);
            Norm             = SubClass_2_2.Compute_ErrorNorms(X.Error,msh);
        end
        
        %% > Step #5.3: Deferred-correction approach.
        % >> Deferred-correction approach: A(H)x=b, where x=? -> A(LO)x=b-[A(HO)x-A(LO)x]_(OLD), where [A(HO)x]_(OLD) and [A(LO)x]_(OLD) are computed explicitly.
        function [X,Norm] = Reconstruct_DC(msh,bnd,blk,F_Vol)
            %  > #1: LO-Compute reconstruction coefficients(Phi,gradPhi).
            [~,~,xf_LO.Phi]     = SubClass_2_2.Reconstruct_Faces(msh,obj.n_LO,1);
            [~,~,xf_LO.gradPhi] = SubClass_2_2.Reconstruct_Faces(msh,obj.n_LO,2);
            xf_LO.Eq            = obj.V.*xf_LO.Phi-obj.Gamma.*xf_LO.gradPhi;
            %  > #2: HO-Compute reconstruction coefficients(Phi,gradPhi).
            [~,~,xf_HO.Phi]     = SubClass_2_2.Reconstruct_Faces(msh,obj.n_HO,1);
            [~,~,xf_HO.gradPhi] = SubClass_2_2.Reconstruct_Faces(msh,obj.n_HO,2);
            xf_HO.Eq            = obj.V.*xf_HO.Phi-obj.Gamma.*xf_HO.gradPhi;
            %  > #3: Assemble matrices(b_LO,a_LO,a_HO,b_HO).
            a_LO                = SubClass_2_2.Assemble_a(msh,obj.n_LO,xf_LO);
            b_LO                = SubClass_2_2.Assemble_b(msh,obj.n_LO,xf_LO,bnd);               
            a_HO                = SubClass_2_2.Assemble_a(msh,obj.n_HO,xf_HO);
            b_HO                = SubClass_2_2.Assemble_b(msh,obj.n_HO,xf_HO,bnd);  
            %  > #4: Assemble RHS.
            X.Phi               = blk.f;
            RHS_1               = b_HO+F_Vol.HO.Ap';
            RHS_2               = a_HO*X.Phi'-b_HO-F_Vol.HO.Ap';
            RHS_3               = a_LO*X.Phi'-b_LO-F_Vol.LO.Ap';
            RHS_T               = RHS_1-(RHS_2-RHS_3);
            %  > #5: Compute approximate solution.
            X.Phi_PDE           = SubClass_2_2.SetUp_bicgstabl(a_LO,RHS_T,10e-12,10e3)';
            %  > #6: Compute error norms.
            X.Error             = abs(X.Phi-X.Phi_PDE);
            Norm                = SubClass_2_2.Compute_ErrorNorms(X.Error,msh);
        end
        
        %% > Step #6: Advance in time.
        
        %% > Step #7: Compute local/global errors.
        function [Norm] = Compute_ErrorNorms(X_Error,msh)
            for i = 1:msh.c.NC
                X(i)            = X_Error(i);
                Norm.E_iX{1}(i) = X(i).*msh.c.Vol(i);
                Norm.E_iX{2}(i) = X(i).^2.*msh.c.Vol(i).^2;
            end
            Norm.E{1} = sum(Norm.E_iX{1})./sum(msh.c.Vol);
            Norm.E{2} = sum(sqrt(Norm.E_iX{2}))./sum(sqrt(msh.c.Vol.^2));
            Norm.E{3} = max(X(i));
        end
        
        %% > Step #8: Compute flux reconstruction SRT.
        function [SRT] = Compute_SRT(msh,bnd,blk,F_Vol)
            if strcmpi(obj.Sim_Type,'Explicit')
                SRT = timeit(@()SubClass_2_2.Reconstruct_ExplicitFlux(msh,bnd,blk,F_Vol));
            elseif strcmpi(obj.Sim_Type,'Implicit')
                SRT = timeit(@()SubClass_2_2.Reconstruct_ImplicitFlux(msh,bnd,blk,F_Vol));
            elseif strcmpi(obj.Sim_Type,'DC')
                SRT = timeit(@()SubClass_2_2.Reconstruct_DC(msh,bnd,blk,F_Vol));
            end
        end
    end
end
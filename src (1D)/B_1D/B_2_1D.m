classdef B_2_1D
    methods (Static)
        %% > Wrap-up B_2 (1D).
        function [msh,pde] = WrapUp_B_2_1D(msh,pde,st,ng,np,v,g,bnd_w,bnd_e)
            % >> Assemble A and B.
            [msh,f,gradf,A,B] = B_2_1D.Assemble_AB(msh,pde,v,g,st,np,ng,bnd_w,bnd_e);
            % >> Compute PDE solution.
            %  > ...cell(s).
            pde = B_2_1D.PDE_cv(msh,pde,A,B);
            %  > ...face(s).
            pde = B_2_1D.PDE_fv(msh,pde,f,gradf,msh.s.bnd.v);
        end
    
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [msh,f,gradf,A,B] = Assemble_AB(msh,pde,v,g,st,np,ng,bnd_w,bnd_e)
            %  > Initialize.
            np_i = ones(1,msh.f.NF);
            np_i = np_i.*np; 
            
            %% > Face 'i'...
            for i = 1:msh.f.NF
                %  > ...stencil indices.
                [msh.s.c{i},msh.s.f{i},msh.s.bnd.i{i}] = ...
                    A_2_1D.SetUp_Stencil_f(msh,i,np_i(i),bnd_w,bnd_e);
                %  > ...stencil coordinates.
                [msh.s.x_v_c{i},msh.s.x_v_f{i},msh.s.x_v_t{i}] = ...
                    A_2_1D.Compute_Coordinates_cft(msh.c.Xc,msh.s.c{i},msh.f.Xv,msh.s.f{i});
                %  > ...Phi_f & GradPhi_f.
                [f{i},gradf{i},xf{i},msh.s.bnd.v{i}] = ...
                    A_2_1D.x_f(pde.sn.f,pde.sn.df,np,msh.s.x_v_t{i},msh.f.Xv(i),msh.s.f{i},msh.s.bnd.i{i},v,g);
            end
            
            %% > Assemble matrix A.
            % >> A (Cell dependent coefficients).
            %  > ... for cell #i: Phi(C) = Phi_f(e)-Phi_f(w).
            %  > Initialize.
            A = zeros(msh.c.NC);
            
            for i = 1:msh.c.NC
                %  > w.
                n_cw       {i}  = length(msh.s.c{i});
                k_cw       {i}  = msh.s.c{i};
                x_fw       {i}  = xf{i};
                A   (i,k_cw{i}) = A(i,k_cw{i})-x_fw{i}(1:n_cw{i});
                %  > e.
                n_ce       {i}  = length(msh.s.c{i+1});
                k_ce       {i}  = msh.s.c{i+1};
                x_fe       {i}  = xf{i+1};
                A   (i,k_ce{i}) = A(i,k_ce{i})+x_fe{i}(1:n_ce{i});
            end
            
            %% > Assemble matrix B.
            % >> B (Face dependent coefficients).
            %  > Initialize.
            B = zeros(msh.c.NC,1);
            
            %  > Boundary contribution(s).
            for i = 1:msh.c.NC
                %  > w.
                if ~isempty(msh.s.f{i})
                    n_fw{i}   = length(msh.s.f{i});
                    k_fw{i}   = n_cw{i}+n_fw{i};
                    vw  {i}   = [msh.s.bnd.v{i}];
                    B   (i,1) = B(i,1)+vw{i}.*x_fw{i}(k_fw{i});
                end
                %  > e.
                if ~isempty(msh.s.f{i+1})
                    n_fe{i}   = length(msh.s.f{i+1});
                    k_fe{i}   = n_ce{i}+n_fe{i};
                    ve  {i}   = [msh.s.bnd.v{i+1}];
                    B   (i,1) = B(i,1)-ve{i}.*x_fe{i}(k_fe{i});
                end
            end
            %  > Source term contribution(s).
            i     = 1:msh.c.NC;
            F_Vol = B_1_2_1D.Compute_SourceTerm(msh,pde,st,ng);
            B     = B+F_Vol(i);
        end
        % >> 1.2. ---------------------------------------------------------
        function [X] = SetUp_Solver(A,B,str)
            if strcmpi(str,'backslash')
                X       = A\B; % Equivalent to: V*((U'*B)./s).
            elseif strcmpi(str,'Tikhonov')
                lmbd    = 1E-06;
                [U,S,V] = svd_lapack(A);
                s       = diag(S);
                X       = V*(s.*(U'*B)./(s.^2+lmbd));
            else
                return;
            end
        end
        % >> 1.3. ---------------------------------------------------------
        %  > 1.3.1. -------------------------------------------------------
        function [pde] = PDE_cv(msh,pde,A,B)
            % >> Compute...
            %  > ...cell nodal solution.
            pde.xn.c = B_2_1D.SetUp_Solver(A,B,'backslash');
            %  > ...cell error/error norms.
            pde.en.c = B_2_1D.Compute_Error(msh,pde,'c');
        end
        %  > 1.3.2. -------------------------------------------------------
        function [pde] = PDE_fv(msh,pde,f,gradf,bnd_v)
            for i = 1:msh.f.NF
                n_c(i) = length(msh.s.c{i});
                k_c{i} = msh.s.c{i};
                v_c{i} = pde.xn.c(k_c{i});
                v_t{i} = v_c{i};
                if ~isempty(msh.s.f{i})
                    n_f(i) = length(msh.s.f{i})+n_c(i);
                    k_f{i} = msh.s.f{1,i};
                    v_f{i} = bnd_v(k_f{i});
                    v_t{i} = [v_t{i};v_f{i}{:}];
                end
                % >> Compute...
                %  > ...face nodal solution.
                pde.xn.f.f (i,1) = f    {i}*v_t{i};
                pde.xn.f.df(i,1) = gradf{i}*v_t{i};
            end
            %  > ...face error.
            pde.en.f = B_2_1D.Compute_Error(msh,pde,'f');
        end
        
        %% > 2. -----------------------------------------------------------
        function [E] = Compute_Error(msh,pde,char_cf)
            switch char_cf
                case 'c'
                    i         = 1:msh.c.NC;
                    Vol (i)   = msh.c.Vol(i);
                    EA  (i)   = abs(pde.sn.c(i)-pde.xn.c(i));
                    E.c (i,1) = EA(i);
                    E1  (i)   = EA(i).*Vol(i);
                    E2  (i)   = E1(i).^2;
                    E.n (1,1) = sum(E1)./sum(Vol);
                    E.n (2,1) = sum(sqrt(E2))./sum(sqrt(Vol));
                    E.n (3,1) = max(EA);
                case 'f'
                    i         = 1:msh.f.NF;
                    E.f (i,1) = abs(pde.sn.f (i)-pde.xn.f.f (i));
                    E.df(i,1) = abs(pde.sn.df(i)-pde.xn.f.df(i));
                otherwise
                    return;
            end
        end
    end
end
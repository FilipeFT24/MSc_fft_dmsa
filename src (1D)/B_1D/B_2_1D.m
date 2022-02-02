classdef B_2_1D
    methods (Static)
        %% > Wrap-up B_2 (1D).
        function [pde] = WrapUp_B_2_1D(msh,pde,ft,st,ng,np,v,g,bnd_w,bnd_e)
            % >> Set...
            %  > ...A and B.
            [A,B] = B_2_1D.Assemble_AB(msh,pde,v,g,st,np,ng,bnd_w,bnd_e);
            %  > ...PDE solution.
            if strcmpi(ft,'Implicit')
                %  > Flux reconstruction: Implicit.
                pde = B_2_1D.PDE_Implicit(msh,pde,A,B);
            elseif strcmpi(ft,'Explicit')
                %  > Flux reconstruction: Explicit.
                pde = B_2_1D.PDE_Implicit(msh,pde,A,B);
            else
                return;
            end
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [A,B] = Assemble_AB(msh,pde,v,g,st,np,ng,bnd_w,bnd_e)
            %% > Face 'i'...
            for i = 1:msh.f.NF
                %  > ...stencil indices.
                [msh.s.c{i},msh.s.f{i}] = ...
                    A_2_1D.SetUp_Stencil_f(i,msh,np);
                %  > ...stencil coordinates.
                [msh.s.x_v_c{i},msh.s.x_v_f{i},msh.s.x_v_t{i}] = ...
                    A_2_1D.Compute_Coordinates_cft(msh,msh.s.c{i},msh.s.f{i});
                %  > ...Phi_f & GradPhi_f.
                [xf.v{i},xf.g{i}] = ...
                    A_2_1D.x_f(v,g,np,msh.s.x_v_t{i},msh.f.Xv(i));
            end
            
            %% > Assemble matrix A.
            % >> A (Cell dependent coefficients).
            %  > ... for cell #i: Phi_C(i) = Phi_f(e)-Phi_f(w).
            %  > Initialize.
            A = zeros(msh.c.NC);
            
            for i = 1:msh.c.NC
                %  > West face.
                n_cw       {i}  = length(msh.s.c{i});
                k_cw       {i}  = msh.s.c{i};
                A   (i,k_cw{i}) = A(i,k_cw{i})-( xf.v{i}  (1:n_cw{i}));
                A   (i,k_cw{i}) = A(i,k_cw{i})-(-xf.g{i}  (1:n_cw{i}));
                %  > East face.
                n_ce       {i}  = length(msh.s.c{i+1});
                k_ce       {i}  = msh.s.c{i+1};
                A   (i,k_ce{i}) = A(i,k_ce{i})+( xf.v{i+1}(1:n_ce{i}));
                A   (i,k_ce{i}) = A(i,k_ce{i})+(-xf.g{i+1}(1:n_ce{i}));
            end
            
            %% > Assemble matrix B.
            % >> B (Face dependent coefficients).
            %  > Initialize.
            B     = zeros(msh.c.NC,1);
            np_2  = np./2;
            %  > Auxiliary array.
            bnd_v = pde.an.bnd;
            
            % >> Boundary contribution(s).
            %  > West boundary.
            k_iw  = 1;
            k_fw  = np_2;
            switch bnd_w
                case 'Dirichlet'
                    for i = k_iw:k_fw
                        B(i,1) = B(i,1)-B_2_1D.Dirichlet_BC(msh.s.f{i},n_cw{i},msh.s.f{i+1},n_ce{i},xf.v{i},xf.g{i},xf.v{i+1},xf.g{i+1},bnd_v);
                    end
                case 'Neumann'
                    for i = k_iw:k_fw
                        B(i,1) = B(i,1)-0;
                    end
                case 'Robin'
                    for i = k_iw:k_fw
                        B(i,1) = B(i,1)-0;
                    end
                otherwise
                    return;
            end
            %  > East boundary.
            k_ie  = msh.c.NC-np_2;
            k_fe  = msh.c.NC;
            switch bnd_e
                case 'Dirichlet'
                    for i = k_ie:k_fe
                        B(i,1) = B(i,1)-B_2_1D.Dirichlet_BC(msh.s.f{i},n_cw{i},msh.s.f{i+1},n_ce{i},xf.v{i},xf.g{i},xf.v{i+1},xf.g{i+1},bnd_v);
                    end
                case 'Neumann'
                    for i = k_ie:k_fe
                        B(i,1) = B(i,1)-0;
                    end
                case 'Robin'
                    for i = k_ie:k_fe
                        B(i,1) = B(i,1)-0;
                    end
                otherwise
                    return;
            end
            % >> Source term contribution(s).
            i     = 1:msh.c.NC;
            F_Vol = B_1_2_1D.Compute_SourceTerm(msh,pde,st,ng);
            B     = B+F_Vol(i);
        end
        %  > 1.1.1. -------------------------------------------------------
        function [Add_t] = Dirichlet_BC(f_w,nc_w,f_e,nc_e,xf_vw,xf_gw,xf_ve,xf_ge,bnd_v)
            %  > Initialize.
            Add_t = 0;
            
            %  > West face.
            if ~isempty(f_w)
                n_fw   =  length(f_w)+nc_w;
                k_fw   =  f_w;
                Add_vw =  xf_vw(n_fw).*bnd_v(k_fw);
                Add_gw = -xf_gw(n_fw).*bnd_v(k_fw);
                Add_tw =  Add_vw+Add_gw;
                Add_t  =  Add_t-Add_tw;
            end
            %  > East face.
            if ~isempty(f_e)
                n_fe   =  length(f_e)+nc_e;
                k_fe   =  f_e;
                Add_ve =  xf_ve(n_fe).*bnd_v(k_fe);
                Add_ge = -xf_ge(n_fe).*bnd_v(k_fe);
                Add_te =  Add_ve+Add_ge;
                Add_t  =  Add_t+Add_te;
            end
        end
        %  > 1.1.2. -------------------------------------------------------
        function [] = Neumann_BC()
        end
         %  > 1.1.3. -------------------------------------------------------
        function [] = Robin_BC()
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
        %  > Implicit flux reconstruction.
        function [pde] = PDE_Implicit(msh,pde,A,B)
            %  > Phi_c.
            pde.Phi = B_2_1D.SetUp_Solver(A,B,'backslash');
            %  > Error norms.
            i       = 1:msh.c.NC;
            X(i)    = abs(pde.an.blk(i)-pde.Phi(i));
            pde.E   = B_2_1D.Compute_ErrorNorms(msh,X);
        end
        % >> 1.4. ---------------------------------------------------------
        %  > Explicit flux reconstruction.
        function [pde] = PDE_Explicit(msh,pde,A,B)
            %  > Error norms.
            X     = abs(A*pde.blk.f-B)';
            pde.E = B_2_2.Compute_ErrorNorms(msh,X);
        end
        
        %% > 2. -----------------------------------------------------------
        function [E]   = Compute_ErrorNorms(msh,X)
            i          = 1:msh.c.NC;
            E.EA(i)    = X(i);
            E_iX{1}(i) = X(i).*msh.c.Vol(i);
            E_iX{2}(i) = X(i).^2.*msh.c.Vol(i).^2;
            E.EN(1)    = sum(E_iX{1})./sum(msh.c.Vol);
            E.EN(2)    = sum(sqrt(E_iX{2}))./sum(sqrt(msh.c.Vol.^2));
            E.EN(3)    = max(X);
        end
    end
end
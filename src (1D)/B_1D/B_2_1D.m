classdef B_2_1D
    methods (Static)
        %% > Wrap-up B_2 (1D).
        function [pde] = WrapUp_B_2_1D(msh,pde,ft,st,ng,np,v,g)
            % >> Set...
            %  > ...A and B.
            [A,B] = B_2_1D.Assemble_AB(msh,pde,v,g,st,np,ng);
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
        function [A,B] = Assemble_AB(msh,pde,v,g,st,np,ng)
            %% > Face 'i'...
            for i = 1:msh.f.NF
                %  > ...stencil indices.
                [msh.s.c{i},msh.s.f{i}] = ...
                    A_2_1D.SetUp_Stencil_f(i,msh,np);
                %  > ...stencil coordinates.
                [msh.s.x_v_c{i},msh.s.x_v_f{i},msh.s.x_v_t{i}] = ...
                    A_2_1D.Compute_Coordinates_cft(msh,msh.s.c{i},msh.s.f{i});
                %  > ...Phi_f & gradPhi_f.
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
            B   = zeros(msh.c.NC,1);
            %  > Auxiliary array.
            bnd = pde.an.bnd;
            
            %  > Compute source term...
            F_Vol = B_1_2_1D.Compute_SourceTerm(msh,pde,st,ng);
            for i = 1:msh.c.NC
                % >> Source term contribution(s).
                B(i,1) = B(i,1)+F_Vol(i);
                % >> Boundary contribution(s).
                %  > West face.
                if ~isempty(msh.s.f{i})
                    n_fw  {i} = length(msh.s.f{i})  +n_cw{i};
                    k_fw  {i} = msh.s.f{i};
                    B   (i,1) = B(i,1)+( xf.v{i}  (n_fw{i})).*bnd(k_fw{i});
                    B   (i,1) = B(i,1)+(-xf.g{i}  (n_fw{i})).*bnd(k_fw{i});
                end
                %  > East face.
                if ~isempty(msh.s.f{i+1})
                    n_fe  {i} = length(msh.s.f{i+1})+n_ce{i};
                    k_fe  {i} = msh.s.f{i+1};
                    B   (i,1) = B(i,1)-( xf.v{i+1}(n_fe{i})).*bnd(k_fe{i});
                    B   (i,1) = B(i,1)-(-xf.g{i+1}(n_fe{i})).*bnd(k_fe{i});
                end
            end
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
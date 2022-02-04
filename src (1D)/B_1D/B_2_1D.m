classdef B_2_1D
    methods (Static)
        %% > Wrap-up B_2 (1D).
        function [msh,pde] = WrapUp_B_2_1D(msh,pde,st,ng,np,v,g,bnd_w,bnd_e)
            % >> Assemble A and B.
            [msh,Tf,A,B] = B_2_1D.Assemble_AB(msh,pde,v,g,st,np,ng,bnd_w,bnd_e);
            % >> Compute PDE solution.
            %  > ...cell(s).
            pde = B_2_1D.PDE_cv(msh,pde,A,B);
            %  > ...face(s).
            pde = B_2_1D.PDE_fv(msh,pde,msh.s.bnd.v,Tf);
        end
    
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [msh,Tf,A,B] = Assemble_AB(msh,pde,v,g,st,np,ng,bnd_w,bnd_e)
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
                %  > ...Tf, xf and boundary values.
                [Tf{i},xf{i},msh.s.bnd.v{i}] = ...
                    A_2_1D.x_f(pde.sn.f(:,1),pde.sn.f(:,2),np_i(i),msh.s.x_v_t{i},msh.f.Xv(i),msh.s.f{i},msh.s.bnd.i{i},v,g);
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
            switch str
                case char('backslash')
                    X       = A\B; % Equivalent to: V*((U'*B)./s).
                case char('Tikhonov')
                    lmbd    = 1E-06;
                    [U,S,V] = svd_lapack(A);
                    s       = diag(S);
                    X       = V*(s.*(U'*B)./(s.^2+lmbd));
                otherwise
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
        function [pde] = PDE_fv(msh,pde,bnd_v,Tf)
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
                j       {i}         = 1:size(Tf{i},1);
                pde.xn.f{i,1}(j{i}) = Tf{i}*v_t{i};
            end
            %  > ...face error.
            pde.en.f = B_2_1D.Compute_Error(msh,pde,'f');
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [E] = Compute_Error(msh,pde,char_cf)
            switch char_cf
                case 'c'
                    %  > Cell(s).
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
                    %  > Fce(s).
                    for i = 1:msh.f.NF
                        sz(i) = length(pde.xn.f{i});
                    end
                    ki = size(pde.fn.f,2);
                    kf = max (sz);
                    if kf > 2
                        [pde.fn.sym(ki+1:kf),pde.fn.f(ki+1:kf)] = ...
                            B_2_1D.Compute_Derivatives_rec(pde.fn.sym{2},ki,kf);
                        i = 1:msh.f.NF;
                        for j = ki+1:kf
                            pde.sn.f(i,j) = pde.fn.f{j}(msh.f.Xv(i));
                        end
                    end
                    for i = 1:msh.f.NF
                        for j = 1:size(pde.xn.f{i},2)
                            E{i,1}(j) = abs(pde.sn.f(i,j)-pde.xn.f{i}(j));
                        end
                    end
                otherwise
                    return;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        function [dnf_sym,dnf_hdl] = Compute_Derivatives_rec(df,ki,kf)
            %  > Symbolic variable.
            syms x;

            for k = 1:kf-ki
                if k == 1
                    dnf_sym{k} = diff(df);
                else
                    dnf_sym{k} = diff(dnf_sym{k-1});
                end
                dnf_hdl{k} = matlabFunction(dnf_sym{k});
            end
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        function [] = Adapt_p()
        end
    end
end
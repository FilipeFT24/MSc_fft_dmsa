classdef B_2_1D
    methods (Static)
        %% > Wrap-up B_2 (1D).
        function [msh,pde] = WrapUp_B_2_1D(msh,pde,st,ng,np,v,g,bnd)
            %% > Initialize...
            %  > ...indices to compute/update stencil,etc.
            stl_p  = repelem(np,msh.f.NF);
            stl_s  = 1:msh.f.NF;
            %  > ...analytic function values (to set bc).
            f(:,1) = pde.sn.f(:,1);
            f(:,2) = pde.sn.f(:,2);
            %  > ...matrices A and B.
            A      = zeros(msh.c.NC);
            B      = zeros(msh.c.NC,1);            
            %% > Set...
            %  > ...source term contribution(s).
            F_Vol = B_1_2_1D.Compute_SourceTerm(msh,pde,st,ng);
            %  > ...matrix B (w/o quadrature).
            i     = 1:msh.c.NC;
            B     = B+F_Vol(i);            
            %% > Update...
            %  > ...stencil and matrices A/B.
            [msh,A,B] = B_2_1D.Update_stl(msh,A,B,stl_p,stl_s,bnd,f,v,g);
            %  > ...cell/face error/error norms.
            pde       = B_2_1D.Update_PDE(msh,pde,A,B,v,g);
        end
    
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [msh,A,B] = Update_stl(msh,A,B,stl_p,stl_s,bnd,f,v,g)
            % >> Compute/update...
            %  > ...stencil.
            msh   = A_2_1D.Problem_SetUp(msh,stl_p,stl_s,bnd,f,v,g);
            xf    = msh.s.xf;
            bnd_s = find(~cellfun(@isempty,msh.s.f));
            bnd_v = cell2mat(msh.s.bnd(bnd_s));
            
            %  > ...matrix A (cell dependent coefficients).
            %  > ...matrix B (face dependent coefficients w/o source term contribution).
            for i = 1:msh.c.NC
                %  > If cell i's west(w) face has been updated...
                if ismembc(i,stl_s)
                    %  > A.
                    kc_w       = msh.s.c{i};
                    xf_w       = xf{i};
                    A(i,kc_w)  = A(i,kc_w)-xf_w(1:length(kc_w));
                    %  > B (boundary contribution(s)).
                    if ismembc(i,bnd_s)
                        vw     = bnd_v(bnd_s == i);
                        B(i,1) = B(i,1)+vw.*xf_w(length(kc_w)+1);
                    end
                end
                %  > If cell i's east(e) face has been updated...
                if ismembc(i+1,stl_s)
                    %  > A.
                    kc_e       = msh.s.c{i+1};
                    xf_e       = xf{i+1};
                    A(i,kc_e)  = A(i,kc_e)+xf_e(1:length(kc_e));
                    %  > B (boundary contribution(s)).
                    if ismembc(i+1,bnd_s)
                        ve     = bnd_v(bnd_s == i+1);
                        B(i,1) = B(i,1)-ve.*xf_e(length(kc_e)+1);
                    end
                end
            end
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
        function [pde] = Update_PDE(msh,pde,A,B,v,g)
            % >> Compute PDE solution...
            %  > ...cell(s).
            pde = B_2_1D.PDE_cv(msh,pde,A,B);
            %  > ...face(s).
            pde = B_2_1D.PDE_fv(msh,pde,v,g);
        end
        %  > 1.3.1. -------------------------------------------------------
        function [pde] = PDE_cv(msh,pde,A,B)
            % >> Compute...
            %  > ...cell nodal solution.
            pde.xn.c = B_2_1D.SetUp_Solver(A,B,'backslash');
            %  > ...cell error/error norms.
            pde.en.c = B_2_1D.Compute_Error(msh,pde,'c');
        end
        %  > 1.3.2. -------------------------------------------------------
        function [pde] = PDE_fv(msh,pde,v,g)
            for i  = 1:msh.f.NF
                kc = msh.s.c{i};
                vc = pde.xn.c(kc);
                vt = vc;
                if ~isempty(msh.s.f{i})
                    vf = msh.s.bnd{i};
                    vt = [vt;vf];
                end
                % >> Compute...
                %  > ...face nodal solution.
                j             = 1:2;
                pde.xn.f(i,j) = vt'*msh.s.Tf{i}(j,:)';
            end
            %  > ...face error/error norms.
            pde.en.f = B_2_1D.Compute_Error(msh,pde,'f',v,g);
        end
        % >> 1.4. ---------------------------------------------------------
        function [E] = Compute_Error(msh,pde,cf,v,g)
            switch cf
                case 'c'
                    % >> Cell(s).
                    i        = 1:msh.c.NC;
                    Vol(i,1) = msh.c.Vol(i);
                    %  > Absolute error.
                    EA (i,1) = abs(pde.sn.c(i,1)-pde.xn.c(i,1));
                    E.c(i,1) = EA(i,1);
                    %  > Error norms.
                    E1 (i,1) = EA(i,1).*Vol(i);
                    E2 (i,1) = E1(i,1).^2;
                    E.n(1,1) = sum(E1)./sum(Vol);
                    E.n(2,1) = sum(sqrt(E2))./sum(sqrt(Vol.^2));
                    E.n(3,1) = max(EA);
                case 'f'
                    % >> Face(s).
                    %  > Absolute error.
                    %  > Column #1: Absolute   error: Phi(f).
                    %  > Column #2: Absolute   error: GradPhi(f).
                    %  > Column #3: Normalized error: [u*(Phi-Phi_A)+g*(GradPhi-GradPhi_A)./Vol]./[u+g./Vol].
                    i        = 1:msh.f.NF;
                    Ls (i,1) = msh.s.Ls;
                    Fvg(i,1) = msh.s.Fvg;
                    E.f(i,1) = abs(pde.sn.f(i,1)-pde.xn.f(i,1));
                    E.f(i,2) = abs(pde.sn.f(i,2)-pde.xn.f(i,2));
                    E.f(i,3) = (v.*E.f(i,1)+g.*E.f(i,2)./Ls(i,1))./Fvg(i,1);
                    %  > Error norms.
                    E1 (i,1) = E.f(i,3).*Ls(i,1);
                    E2 (i,1) = E1(i,1).^2;
                    E.n(1,1) = sum(E1)./sum(Ls);
                    E.n(2,1) = sum(sqrt(E2))./sum(sqrt(Ls.^2));
                    E.n(3,1) = max(E.f(:,3));
                otherwise
                    return;
            end
        end
    end
end
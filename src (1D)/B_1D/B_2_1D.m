classdef B_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [s,A,B] = Update_stl(msh,stl_p,stl_s,A,B,sn,bnd,v,g)
            % >> Compute/update...
            %  > ...stencil.
            s     = A_2_1D.Problem_SetUp(msh,stl_p,stl_s,sn,bnd,v,g);
            xf    = s.xf;
            bnd_s = find(~cellfun(@isempty,s.f));
            bnd_v = cell2mat(s.bnd(bnd_s));
            
            %  > ...matrix A (cell dependent coefficients).
            %  > ...matrix B (face dependent coefficients w/o source term contribution).
            for i = 1:msh.c.NC
                %  > If cell i's west(w) face has been updated...
                if ismembc(i,stl_s)
                    %  > A.
                    kc_w       = s.c{i};
                    xf_w       = xf {i};
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
                    kc_e       = s.c{i+1};
                    xf_e       = xf {i+1};
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
        function [en] = Update_pde(msh,A,B,sn,s,v,g)
            % >> Update pde...
            %  > ...cell(s).
            xn.c = A\B;
            %  > ...face(s).
            for i  = 1:msh.f.NF
                kc = s.c{i};
                vc = xn.c(kc);
                vt = vc;
                if ~isempty(s.f{i})
                    vf = s.bnd{i};
                    vt = [vt;vf];
                end
                j         = 1:2;
                xn.f(i,j) = vt'*s.Tf{i}(j,:)';
            end
            %  > ...cell/face error/error norms.
            [en.c,en.f] = B_2_1D.Compute_Error(msh,sn,xn,s,v,g);
        end
        % >> 1.3. ---------------------------------------------------------
        function [Ec,Ef] = Compute_Error(msh,sn,xn,s,v,g)
            % >> Cell(s).
            i         = 1:msh.c.NC;
            Vol(i,1)  = msh.c.Vol(i);
            %  > Absolute error.
            Ec.c(i,1) = abs(sn.c(i,1)-xn.c(i,1));
            %  > Error norms.
            Ec_1(i,1) = Ec.c(i,1).*Vol(i);
            Ec_2(i,1) = Ec_1(i,1).^2;
            Ec.n(1,1) = sum(Ec_1)./sum(Vol);
            Ec.n(2,1) = sum(sqrt(Ec_2))./sum(sqrt(Vol.^2));
            Ec.n(3,1) = max(Ec_1);
            
            % >> Face(s).
            i         = 1:msh.f.NF;
            Ls (i,1)  = s.Ls;
            Fvg(i,1)  = s.Fvg;
            %  > Absolute error: 1) Column #1: Absolute   error: Phi(f).
            %                    2) Column #2: Absolute   error: GradPhi(f).
            %                    3) Column #3: Normalized error: [u*(Phi-Phi_A)+g*(GradPhi-GradPhi_A)./Vol]./[u+g./Vol].
            Ef.f(i,1) = abs(sn.f(i,1)-xn.f(i,1));
            Ef.f(i,2) = abs(sn.f(i,2)-xn.f(i,2));
            Ef.f(i,3) = (v.*Ef.f(i,1)+g.*Ef.f(i,2)./Ls(i,1))./Fvg(i,1);
            %  > Error norms.
            Ef_1(i,1) = Ef.f(i,3).*Ls(i,1);
            Ef_2(i,1) = Ef_1(i,1).^2;
            Ef.n(1,1) = sum(Ef_1)./sum(Ls);
            Ef.n(2,1) = sum(sqrt(Ef_2))./sum(sqrt(Ls.^2));
            Ef.n(3,1) = max(Ef.f(:,3));
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [msh,pde] = SetUp_Problem(msh,pde,v,g,bnd,np,p_adapt)
            %  > Initialize.
            A  = zeros(msh.c.NC);
            B  = zeros(msh.c.NC,1);
            
            %  > Select...
            switch p_adapt
                case false
                    %  > Initialize/add...
                    %  > ...indices to compute/update stencil,etc.
                    stl_p   = repelem(np,msh.f.NF);
                    stl_s   = 1:msh.f.NF;
                    %  > ...source term contribution.
                    B = B+pde.FV;
                    
                    %  > Solve PDE.
                    [s,A,B] = B_2_1D.Update_stl(msh,stl_p,stl_s,A,B,pde.sn,bnd,v,g);
                    [e]     = B_2_1D.Update_pde(msh,A,B,pde.sn,s,v,g);
                    %  > Update structure.
                    msh.s   = s;
                    pde.e   = e;
                case true
                    
                    
                    s=1;
                    
                    
                otherwise
                    return;
            end
        end
    end
end
classdef B_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [s,xn] = Update_stl(msh,s,stl_p,stl_s,stl_t,A,B,sn,bnd,v,g)
            % >> Compute/update...
            %  > ...stencil.
            s  = A_2_1D.Stencil_SetUp(msh,s,stl_p,stl_s,stl_t,sn,bnd,v,g);
            xf = s.xf;
            for i = 1:size(s.f,1)
                bnd_s{i} = find(~cellfun(@isempty,s.f(i,:)));
                bnd_v{i} = [s.bnd{i,bnd_s{i}}];   
            end
            %  > ...sign (convection/diffusion)
            sign_x = [1,-1];
            
            %  > ...matrix A (cell dependent coefficients).
            %  > ...matrix B (face dependent coefficients w/o source term contribution).
            for i = 1:msh.c.NC
                for j = 1:size(s.xf,1)
                    % >> West face contribution.
                    %  > A.
                    kc_w    {j}  = s.c{j,i};
                    xf_w    {j}  = xf {j,i};
                    A(i,kc_w{j}) = A(i,kc_w{j})-sign_x(j).*xf_w{j}(1:length(kc_w{j}));
                    %  > B.
                    if ismembc(i,bnd_s{j})
                        fw  {j}  = bnd_v{j}(bnd_s{j} == i);
                        B (i,1)  = B(i,1)+sign_x(j).*fw{j}.*xf_w{j}(length(kc_w{j})+1);
                    end
                    % >> East face contribution.
                    %  > A.
                    kc_e    {j}  = s.c{j,i+1};
                    xf_e    {j}  = xf {j,i+1};
                    A(i,kc_e{j}) = A(i,kc_e{j})+sign_x(j).*xf_e{j}(1:length(kc_e{j}));
                    %  > B.
                    if ismembc(i+1,bnd_s{j})
                        fe  {j}  = bnd_v{j}(bnd_s{j} == i+1);
                        B (i,1)  = B(i,1)-sign_x(j).*fe{j}.*xf_e{j}(length(kc_e{j})+1);
                    end
                end
            end
            xn.c = A\B;
        end
        % >> 1.2. ---------------------------------------------------------
        function [e,x] = Update_pde(msh,x,sn,s)
            % >> Update pde...
            %  > ...face(s).
            for i  = 1:msh.f.NF
                for j = 1:size(s.xf,1)
                    if ~isempty(s.c{j,i})
                        kc{j} = s.c{j,i};
                        vc{j} = x.c(kc{j});
                        vt{j} = vc {j};
                        if ~isempty(s.f{j,i})
                            vf{j} = s.bnd{j,i};
                            vt{j} = [vt{j};vf{j}];
                        end
                    else
                        vt{j} = s.bnd{j,i};
                    end
                    x.f(i,j) = s.Tf{j,i}*vt{j};
                end
            end
            %  > ...cell/face error/error norm(s).
            [e.c,e.f] = B_2_1D.Compute_Error(msh,sn,x,s);
        end
        % >> 1.3. ---------------------------------------------------------
        function [Ec,Ef] = Compute_Error(msh,sn,xn,s)
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
            i        = 1:msh.f.NF;
            for j = 1:size(xn.f,2)
                %  > Stencil width.
                Ls (i,j) = s.Ls(j,i);
                %  > Absolute error: 1) Column #1: Absolute   error: Phi(f).
                %                    2) Column #2: Absolute   error: GradPhi(f).
                Ef.f(i,j) = abs(sn.f(i,j)-xn.f(i,j));
                %  > Error norms.
                Ef_1(i,j) = Ef.f(i,j).*Ls(i,j);
                Ef_2(i,j) = Ef_1(:,j).^2;
                Ef.n(1,j) = sum(Ef_1(:,j))./sum(Ls(:,j));
                Ef.n(2,j) = sum(sqrt(Ef_2(:,j)))./sum(sqrt(Ls(:,j).^2));
                Ef.n(3,j) = max(Ef_1(:,j));
            end
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [stl] = Select_f(stl,e)
            %  > Size.
            [m,n] = size(e.f.f);
            
            %  > Check rules...
            for j = 1:n
                k = 1;
                for i = 1:m
                    %  > Increase polynomial order by +1.
                    if e.f.f(i,j) > e.f.n(1,j)
                        stl.p(j,i) = stl.p(j,i)+1;
                        stl_s(j,k) = i;
                        k          = k+1;
                        stl.t(j,i) = "UDS";
                    end
                end
                %  > Check "Neighbour rule" (p-version).
                for i = 1:m
                    if (i > 1 && i < m) && (stl.p(j,i) < stl.p(j,i-1) && stl.p(j,i) < stl.p(j,i+1))
                        stl.p(j,i) = stl.p(j,i)+1;
                        stl_s(j,k) = i;
                        k          = k+1;
                        stl.t(j,i) = "UDS";
                    end
                end
                stl.s = sort(stl_s,2);
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        function [] = Rule_1()
        end
        %  > 2.2.2. -------------------------------------------------------
        function [] = Rule_2()
        end 
    end
end
classdef B_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [s,x,ea] = Update_stl(msh,s,stl_p,stl_s,stl_t,A,B,a,bnd,v,g)
            % >> Compute/update...
            %  > ...stencil.
            s  = A_2_1D.Stencil_SetUp(msh,s,stl_p,stl_s,stl_t,a,bnd,v,g);
            %  > ...xf.
            vg = [v,g];
            for j = 1:size(s.f,1)
                for i = 1:size(s.f,2)
                    xf{j,i} = vg(j).*s.xf{j,i};
                end
                bnd_s{j} = find(~cellfun(@isempty,s.f(j,:)));
                bnd_v{j} = [s.bnd{j,bnd_s{j}}];   
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
            x.c = A\B;
            ea  = abs(A*a.c(:,1)-B);
        end
        % >> 1.2. ---------------------------------------------------------
        function [e,x] = Update_pde(msh,a,s,x,ea,v,g)
            % >> Update pde...
            %  > ...face(s).
            for i  = 1:msh.f.NF
                for j = 1:size(s.xf,1)
                    if ~isempty(s.c{j,i})
                        kc{j} = s.c{j,i};
                        ca{j} = a.c(kc{j},1);
                        va{j} = ca {j};
                        cx{j} = x.c(kc{j},1);
                        vx{j} = cx {j};
                        if ~isempty(s.f{j,i})
                            vf{j} = s.bnd{j,i};
                            va{j} = [va{j};vf{j}];
                            vx{j} = [vx{j};vf{j}];
                        end
                    else
                        va{j}  = s.bnd{j,i};
                        vx{j}  = s.bnd{j,i};
                    end
                    x.f.a(i,j) = s.xf{j,i}*va{j};
                    x.f.f(i,j) = s.xf{j,i}*vx{j};
                end
            end
            %  > ...cell/face error/error norm(s).
            [e.c,e.f] = B_2_1D.Compute_Error(msh,a,x,s,ea,v,g);
        end
        % >> 1.3. ---------------------------------------------------------
        function [Ec,Ef] = Compute_Error(msh,a,x,s,ea,v,g)
            % >> Cell(s).
            vg        = [v,g];
            i         = 1:msh.c.NC;
            Vol(i,1)  = msh.c.Vol(i);
            %  > Absolute error.
            Ec.c.a(i,1) = ea;
            Ec.c.c(i,1) = abs(a.c(i,1)-x.c(i,1));
            %  > Error norms.
            Ec_1.a(i,1) = Ec.c.a(i,1).*Vol(i);
            Ec_1.c(i,1) = Ec.c.c(i,1).*Vol(i);
            Ec_2.a(i,1) = Ec_1.a(i,1).^2;
            Ec_2.c(i,1) = Ec_1.c(i,1).^2;
            Ec.n.a(1,1) = sum(Ec_1.a)./sum(Vol);
            Ec.n.a(2,1) = sum(sqrt(Ec_2.a))./sum(sqrt(Vol.^2));
            Ec.n.a(3,1) = max(Ec_1.a);
            Ec.n.c(2,1) = sum(sqrt(Ec_2.c))./sum(sqrt(Vol.^2));
            Ec.n.c(1,1) = sum(Ec_1.c)./sum(Vol);
            Ec.n.c(3,1) = max(Ec_1.c);
            
            % >> Face(s).
            %  > Absolute error: 1. Column #1: Absolute error: phi_f.
            %  >                 2. Column #2: Absolute error: gradphi_f.
            i = 1:msh.f.NF;
            for j = 1:size(x.f.a,2)
                %  > Absolute error.
                Ls    (i,j) = s.Ls(j,i);
                Ef.f.a(i,j) = abs((a.f(i,j)-x.f.a(i,j))).*vg(j);
                Ef.f.f(i,j) = abs(a.f(i,j)-x.f.f(i,j));
                %  > Error norms.
                Ef_1.a(i,j) = Ef.f.a(i,j).*Ls(i,j);
                Ef_1.f(i,j) = Ef.f.f(i,j).*Ls(i,j);
                Ef_2.a(i,j) = Ef_1.a(:,j).^2;
                Ef_2.f(i,j) = Ef_1.f(:,j).^2;
                Ef.n.a(1,j) = sum(Ef_1.a(:,j))./sum(Ls(:,j));
                Ef.n.a(2,j) = sum(sqrt(Ef_2.a(:,j)))./sum(sqrt(Ls(:,j).^2));
                Ef.n.a(3,j) = max(Ef_1.a(:,j)); 
                Ef.n.f(1,j) = sum(Ef_1.f(:,j))./sum(Ls(:,j));
                Ef.n.f(2,j) = sum(sqrt(Ef_2.f(:,j)))./sum(sqrt(Ls(:,j).^2));
                Ef.n.f(3,j) = max(Ef_1.f(:,j)); 
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Switch_EE(ee)
            switch ee
                case "1"
                    %  > Select low-order(LO) and high-order(HO) schemes.
                    LO  = 2;
                    HO  = 4;
                    %  > Compute error estimation...
                    pde = B_2_1D.EE_1(msh,v,g,ft,LO,HO);
                otherwise
                    return;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        function [] = EE_1(msh,v,g,ft,LO,HO)
            [pde,A,B,stl,s] = B_1D.Initialize(msh,v,g,ft,t1,t2);
        end
    end
end
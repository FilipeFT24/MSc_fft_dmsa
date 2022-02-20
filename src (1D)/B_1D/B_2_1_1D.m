classdef B_2_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [stl,s,x,et_c] = Update_stl(msh,stl,s,A,B,a,bnd,v,g)
            % >> Compute/update...
            %  > ...stencil.
            [stl,s] = A_2_1D.SetUp_stl(msh,stl,s,a,bnd,v,g);
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
            x.c  = A\B;
            et_c = A*a.c(:,1)-B;
        end
        % >> 1.2. ---------------------------------------------------------
        function [e,x] = Update_pde(msh,a,s,x,et_c,v,g)
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
            [e.c,e.f,e.t] = B_2_1_1D.Compute_Error(msh,a,x,et_c,v,g);
        end
        % >> 1.3. ---------------------------------------------------------
        function [Ec,Ef,Et] = Compute_Error(msh,a,x,et_c,v,g)
            %  > Auxiliary variables.
            [m,n] = size(x.f.f);
            vg    = [v,g];
            sign  = [-1,1];

            % >> Face(s).
            for j = 1:n
                for i = 1:m
                    %  > Absolute error.
                    ef_f(i,j) = a.f(i,j)-x.f.f(i,j);
                    ef_t(i,j) = (a.f(i,j)-x.f.a(i,j)).*vg(j);
                    Ef.f(i,j) = abs(ef_f(i,j));
                    Et.f(i,j) = abs(ef_t(i,j));
                end
                Ef.n  (1,j) = mean(Ef.f(:,j));
                Et.n.f(1,j) = mean(Et.f(:,j));
            end
            Et.n.t(1,1) = 1./2.*(Et.n.f(1,1)+Et.n.f(1,2));
            
            % >> Cell(s).
            for i = 1:m-1
                %  > Absolute error.
                ec_c  (i,1) = a.c(i,1)-x.c(i,1);
                Ec.c  (i,1) = abs(ec_c(i));
                %  > Error transported by convection/diffusion.
                for j = 1:n+1
                    if j == 1
                        Et.c(i,j) = abs(et_c(i));
                    else
                        Et.c(i,j) = abs(ef_t(i+1,j-1)-ef_t(i,j-1));
                    end
                end
                Vol   (i,1) = msh.c.Vol(i);
                %  > Error norms.
                Ec_1.c(i,1) = Ec.c  (i).*Vol(i);
                Ec_1.t(i,1) = Et.c  (i).*Vol(i);
                Ec_2.c(i,1) = Ec_1.c(i).^2;
                Ec_2.t(i,1) = Ec_1.t(i).^2;
            end
            Ec.n  (1,1) = sum(Ec_1.c)./sum(Vol);
            Ec.n  (2,1) = sum(sqrt(Ec_2.c))./sum(sqrt(Vol.^2));
            Ec.n  (3,1) = max(Ec.c);
            Et.n.c(1,1) = sum(Ec_1.t)./sum(Vol);
            Et.n.c(2,1) = sum(sqrt(Ec_2.t))./sum(sqrt(Vol.^2));
            Et.n.c(3,1) = max(Et.c(:,1));
        end
    end
end
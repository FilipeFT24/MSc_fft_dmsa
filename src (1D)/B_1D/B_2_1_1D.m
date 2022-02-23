classdef B_2_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [stl,s,A,B,x,et_c] = Update_stl(msh,stl,s,A,B,a,st,bnd,v,g)
            %  > Auxiliary variables.
            NC = msh.c.NC;
            k  = 1:size(stl.s,2);
            vg = [v,-g];
            Ac = zeros(NC);
            Bc = zeros(NC,1);
                                   
            % >> Compute/update...
            %  > ...stencil.
            [stl,s] = A_2_1D.SetUp_stl(msh,stl,s,a,bnd,v,g);
            %  > ...xf.        
            for i = 1:size(s.f,1)
                bnd_s{i} = find(~cellfun(@isempty,s.f(i,:)));
                bnd_v{i} = [s.bnd{i,bnd_s{i}}];   
            end
            %  > ...A/B.
            for i = k
                for j = 1:NC
                    %  > Initialize.
                    [u,A{i}(j,:),B{i}(j,1)] = ...
                        B_2_1_1D.Initialize(j,stl.s{i},s.c(i,:),s.xf(i,:),bnd_s{i},bnd_v{i},NC,A{i}(j,:),B{i}(j,1));
                    
                    %  > Update A.
                    if u.c.uA
                        A{i}(j,u.w.kc) = B_2_1_1D.Assemble_A("w",A{i}(j,u.w.kc),u.w.Af);
                        A{i}(j,u.e.kc) = B_2_1_1D.Assemble_A("e",A{i}(j,u.e.kc),u.e.Af); 
                    end
                    %  > Update B.
                    if u.w.uB
                        B{i}(j,1) = B_2_1_1D.Assemble_B("w",B{i}(j,1),u.w.fv,u.w.Bf);
                    end
                    if u.e.uB
                        B{i}(j,1) = B_2_1_1D.Assemble_B("e",B{i}(j,1),u.e.fv,u.e.Bf);
                    end
                end
                %  > Add convective/diffusive contributions (cumulative matrices)...
                Ac = Ac+vg(i).*A{i};
                Bc = Bc+vg(i).*B{i};
            end
            %  > Add source term...
            Bc = Bc+st;
                       
            %  > Update solution(s)...
            x.c  = Ac\Bc;
            et_c = Ac*a.c(:,1)-Bc;
        end
        %  > 1.1.1. -------------------------------------------------------
        %  > Initialize matrix coefficients.
        function [u,A,B] = Initialize(j,stl_s,sc,xf,bnd_s,bnd_v,NC,A,B)
            %  > Initialize...
            u.c.uA = false;
            u.w.uB = false;
            u.e.uB = false;
            
            if any(ismembc([j,j+1],stl_s))
                % >> A.
                %  > Update...
                A      = zeros(1,NC);
                u.c.uA = true;
                %  > West face contribution.
                u.w.kc = sc{j};
                l.w    = length(u.w.kc);
                u.w.Af = xf{j}(1:l.w);
                %  > East face contribution.
                u.e.kc = sc{j+1};
                l.e    = length(u.e.kc);
                u.e.Af = xf{j+1}(1:l.e);

                % >> B.
                if any(ismembc([j,j+1],bnd_s))
                    %  > Update...
                    B  = 0;
                    if ismembc(j,bnd_s)
                        u.w.uB = true;
                        u.w.fv = bnd_v(bnd_s == j);
                        u.w.Bf = xf{j}(l.w+1);
                    end
                    if ismembc(j+1,bnd_s)
                        u.e.uB = true;
                        u.e.fv = bnd_v(bnd_s == j+1);
                        u.e.Bf = xf{j+1}(l.e+1);
                    end
                end
            end
        end
        %  > 1.1.2. -------------------------------------------------------
        %  > Update matrix A (cell dependent coefficients).
        function [A] = Assemble_A(str_f,A,xf)
            switch str_f
                case "w"
                    %  > West face contribution.
                    A = A-xf;
                case "e"
                    %  > East face contribution.
                    A = A+xf;
                otherwise
                    return;
            end
        end
        %  > 1.1.3. -------------------------------------------------------
        %  > Update matrix B (face dependent coefficients w/o source term contribution).
        function [B] = Assemble_B(str_f,B,f,xf)
            switch str_f
                case "w"
                    %  > West face contribution.
                    B = B+f.*xf;
                case "e"
                    %  > East face contribution.
                    B = B-f.*xf;
                otherwise
                    return;
            end
        end
        % >> 1.2. ---------------------------------------------------------
        function [e,x] = Update_pde(msh,a,s,x,et_c,v,g)
            % >> Update pde...
            %  > ...face(s).
            for i  = 1:msh.f.NF
                for j = 1:size(s.xf,1)
                    if ~isempty(s.c{j,i})
                        %  > Cell indices.
                        kc{j} = s.c{j,i};
                        %  > Analytic/PDE values.
                        va{j} = a.c(kc{j},1);
                        vx{j} = x.c(kc{j},1);
                        
                        %  > Add boundary contribution?...
                        if ~isempty(s.f{j,i})
                            vf{j} = s.bnd{j,i};
                            va{j} = [va{j};vf{j}];
                            vx{j} = [vx{j};vf{j}];
                        end
                    else
                        va{j} = s.bnd{j,i};
                        vx{j} = s.bnd{j,i};
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
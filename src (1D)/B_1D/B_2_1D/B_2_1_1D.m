classdef B_2_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [pde,s,stl] = Update_pde(obj,msh,pde,s,stl)
            %  > Auxiliary variables.
            NC = msh.c.NC;
            k  = 1:size(stl.s,2);
            vg = [obj.v,-obj.g];
            Ac = zeros(NC);
            Bc = zeros(NC,1);
            
            % >> Compute/update...
            %  > ...stencil.
            [stl,s] = A_2_1D.SetUp_stl(obj,msh,pde,s,stl);
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
                        B_2_1_1D.Initialize(j,NC,stl.s{i},s.c(i,:),s.xf(i,:),s.A{i}(j,:),s.B{i}(j,1),bnd_s{i},bnd_v{i});
                    
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
                %  > Add convective/diffusive contributions (cumulative matrices).
                Ac = Ac+vg(i).*A{i};
                Bc = Bc+vg(i).*B{i};
            end
            %  > Add source term.
            Bc  = Bc+pde.f.st;
            
            %  > Update problem variables.
            s.A = A;
            s.B = B;
            pde = B_2_1_1D.Compute_Error(obj,msh,pde,s,Ac,Bc);
        end
        %  > 1.1.1. -------------------------------------------------------
        %  > Initialize matrix coefficients.
        function [u,A,B] = Initialize(j,NC,stl_s,sc,xf,A,B,bnd_s,bnd_v)
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
        %  > Compute error/error norms.
        function [pde] = Compute_Error(obj,msh,pde,s,A,B)
            %  > Auxiliary variables.
            vg = [obj.v,obj.g];
            l  = size(s.xf,1);
            m  = msh.f.NF;
            n  = length(vg);
            
            %  > Nodal solution.
            x.c = A\B;
            %  > Face (reconstructed) values.
            for i = 1:m
                for j = 1:l
                    if ~isempty(s.c{j,i})
                        %  > Analytic/PDE values.
                        kc{j} = s.c{j,i};
                        va{j} = pde.a.c(kc{j},1);
                        vx{j} = x.c    (kc{j},1);
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
                    x.v.va{j,i} = va{j};
                    x.v.vx{j,i} = vx{j};
                    x.f.a (i,j) = s.xf{j,i}*va{j};
                    x.f.f (i,j) = s.xf{j,i}*vx{j};
                end
            end
            
            % >> Face/cell absolute error/truncation error.
            %  > Face(s).
            for j = 1:n
                i = 1:m;
                %  > Error distribution.
                ef_f  (i,j) = pde.a.f(i,j)-x.f.f(i,j);
                et_f  (i,j) = vg(j).*(pde.a.f(i,j)-x.f.a(i,j));
                %  > Absolute/truncation error distribution.
                ef.f  (i,j) = abs(ef_f(i,j));
                et.f  (i,j) = abs(et_f(i,j));
                %  > Absolute/truncation error mean norm.
                ef.n  (1,j) = mean(ef.f(i,j));
                et.n.f(1,j) = mean(et.f(i,j));
            end
            et.n.t(1,1) = 1./2.*(et.n.f(1,1)+et.n.f(1,2));
            
            %  > Cell(s).
            et_c = A*pde.a.c(:,1)-B;
            for i = 1:m-1
                %  > Absolute error distribution.
                ec.c  (i,1) = abs(pde.a.c(i,1)-x.c(i,1));
                %  > Absolute error transported by convection/diffusion.
                for j = 1:n+1
                    if j == 1
                        et.c(i,j) = abs(et_c(i));
                    else
                        et.c(i,j) = abs(et_f(i+1,j-1)-et_f(i,j-1));
                    end
                end
                %  > Cell volume.
                Vol   (i,1) = msh.c.Vol(i);
                %  > Absolute/truncation error norms.
                ec_1.c(i,1) = ec.c  (i).*Vol(i);
                ec_1.t(i,1) = et.c  (i).*Vol(i);
                ec_2.c(i,1) = ec_1.c(i).^2;
                ec_2.t(i,1) = ec_1.t(i).^2;
            end
            ec.n  (1,1) = sum(ec_1.c)./sum(Vol);
            ec.n  (2,1) = sum(sqrt(ec_2.c))./sum(sqrt(Vol.^2));
            ec.n  (3,1) = max(ec.c);
            et.n.c(1,1) = sum(ec_1.t)./sum(Vol);
            et.n.c(2,1) = sum(sqrt(ec_2.t))./sum(sqrt(Vol.^2));
            et.n.c(3,1) = max(et.c(:,1));
            
            %  > Update 'pde'...
            pde.x   = x;
            pde.e.c = ec;
            pde.e.f = ef;
            pde.e.t = et;
        end
    end
end
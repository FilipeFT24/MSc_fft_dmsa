classdef B_2_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [pde,s,stl] = WrapUp_B_2_1_1D(obj,msh,pde,s,stl)
            %  > Update 's' and 'stl' structures and create 'x' field.
            [s,stl,x.c]   = B_2_1_1D.Update_s(obj,msh,pde,s,stl);
            %  > Update 'x' and 'e' fields.
            [x.f.a,x.v.a] = B_2_1_1D.Update_pde_x (msh,s,pde.a);
            [x.f.f,x.v.f] = B_2_1_1D.Update_pde_x (msh,s,x);
            [e.c]         = B_2_1_1D.Update_pde_ec(pde.a,x,msh.c.Vol);
            [e.f]         = B_2_1_1D.Update_pde_ef(pde.a,x);
            [e.t]         = B_2_1_1D.Update_pde_et(obj,msh,s,pde.a,x); 
            %  > Set 'pde' fields.
            pde.x         = x;
            pde.e         = e;
        end
               
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
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
        % >> 2.2. ---------------------------------------------------------
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
        % >> 2.3. ---------------------------------------------------------
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
        % >> 2.4. ---------------------------------------------------------
        %  > Update 's' structure.
        function [s,stl,nc] = Update_s(obj,msh,pde,s,stl)
            %  > Auxiliary variables.
            NC = msh.c.NC;
            k  = 1:size(stl.s,2);
            vg = [obj.v,-obj.g];
            Ac = zeros(NC);
            Bc = zeros(NC,1);
            
            %  > Update stencil.
            [s,stl] = A_2_1D.SetUp_stl(obj,msh,pde,s,stl);
            %  > Update  xf.
            for i = 1:size(s.f,1)
                bnd_s{i} = find(~cellfun(@isempty,s.f(i,:)));
                bnd_v{i} = [s.bnd{i,bnd_s{i}}];
            end
            %  > Update A/B.
            for i = k
                for j = 1:NC
                    %  > Initialize.
                    [u,A{i}(j,:),B{i}(j,1)] = ...
                        B_2_1_1D.Initialize(j,NC,stl.s{i},s.c(i,:),s.xf(i,:),s.A{i}(j,:),s.B{i}(j,1),bnd_s{i},bnd_v{i});
                    
                    %  > A.
                    if u.c.uA
                        A{i}(j,u.w.kc) = B_2_1_1D.Assemble_A("w",A{i}(j,u.w.kc),u.w.Af);
                        A{i}(j,u.e.kc) = B_2_1_1D.Assemble_A("e",A{i}(j,u.e.kc),u.e.Af);
                    end
                    %  > B.
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
            Bc   = Bc+pde.f.st;
            %  > Nodal solution.
            nc   = Ac\Bc;
            %  > Update fields.
            s.A  = A;
            s.B  = B;
            s.Ac = Ac;
            s.Bc = Bc;
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > Update analytic 'pde.x.(...)' field.
        function [rf,vf] = Update_pde_x(msh,s,n)
            %  > Auxiliary variables.
            l = size(s.xf,1);
            m = msh.f.NF;
            
            %  > Face (reconstructed) values.
            for i = 1:m
                for j = 1:l
                    if ~isempty(s.c{j,i})
                        %  > Nodal indices/values.
                        kc{j} = s.c{j,i};
                        vc{j} = n.c(kc{j},1);
                        %  > Add boundary contribution?...
                        if ~isempty(s.f{j,i})
                            vb{j} = s.bnd{j,i};
                            vf{j} = [vc{j};vb{j}];
                        end
                    else
                        vf{j} = s.bnd{j,i};
                    end
                    rf(i,j) = s.xf{j,i}*vf{j};
                end
            end
        end
        % >> 3.2. ---------------------------------------------------------
        %  > 3.2.1. -------------------------------------------------------
        %  > Update 'pde.e.c' field (cell error).
        function [ec] = Update_pde_ec(a,x,Vol)
            %  > Absolute error distribution.
            ec.c(:,1) = abs(a.c(:,1)-x.c(:,1));
            %  > Error norms.
            ec_1(:,1) = ec.c.*Vol;
            ec_2(:,1) = ec_1.^2;
            ec.n(1,1) = sum(ec_1)./sum(Vol);
            ec.n(2,1) = sum(sqrt(ec_2))./sum(sqrt(Vol.^2));
            ec.n(3,1) = max(ec.c);
        end
        %  > 3.2.2. -------------------------------------------------------
        %  > Update 'pde.e.f' field (face error).
        function [ef] = Update_pde_ef(a,x)
            %  > Convective/diffusive components.
            i = 1:2;
            %  > Absolute error distribution.
            ef.f(:,i) = abs(a.f(:,i)-x.f.f(:,i));
            %  > Mean (absolute) error.
            ef.n(1,i) = mean(ef.f(:,i));
        end
        %  > 3.2.3. -------------------------------------------------------
        %  > Update 'pde.e.t' field (truncation error).
        function [et] = Update_pde_et(obj,msh,s,a,x)
            %  > Weighted (convective/diffusive) components.
            vg = [obj.v,obj.g];
            nc = length(vg);
            i  = 1:nc;
            j  = 1:nc+1;
            
            % >> Face(s).
            %  > Truncation error distribution.
            et_f  (:,i) = a.f(:,i)-x.f.a(:,i);
            %  > Weighted truncation error distribution.
            et_f  (:,i) = vg(i).*et_f(:,i);
            %  > Absolute weighted truncation error distribution.
            et.f  (:,i) = abs(et_f(:,i));
            %  > Mean (absolute weighted) truncation error.
            et.n.f(:,i) = mean(et.f(:,i));
            %  > Reference mean (absolute weighted) truncation error.
            et.n.t      = mean(et.n.f);
            
            % >> Cell(s).
            %  > Truncation error distribution.
            et_c  = s.Ac*a.c(:,1)-s.Bc;
            %  > Absolute truncation error transported by convection and/or diffusion.
            for k = j
                if k == 1
                    %  > w/ convection and diffusion.
                    et.c(:,k) = abs(et_c);
                else
                    for l = 1:msh.c.NC
                        %  > w/ convection or diffusion.
                        et.c(i,k) = abs(et_f(l+1,k-1)-et_f(l,k-1));
                    end
                end
            end
            %  > Error norms.
            Vol         = msh.c.Vol;
            ec_1        = et.c.*Vol;
            ec_2        = ec_1.^2;
            et.n.c(1,:) = sum(ec_1,1)./sum(Vol);
            et.n.c(2,:) = sum(sqrt(ec_2),1)./sum(sqrt(Vol.^2));
            et.n.c(3,:) = max(et.c);
        end
    end
end
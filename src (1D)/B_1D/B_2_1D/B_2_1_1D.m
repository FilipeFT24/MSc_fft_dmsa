classdef B_2_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [msh,pde] = SetUp_p(inp,msh,pde,s,stl)
            %  > Update 's' structure and 'x' field.
            s     = B_2_1_1D.Update_s     (inp,msh,pde,s,stl);
            x.c   = B_2_1_1D.Update_xc    (s);
            %  > Compute cell/face errors.
            e.f   = B_2_1_1D.Update_ef    (s,x.c,pde.av.f);
            e.c   = B_2_1_1D.Update_ec    (x.c,pde.av.c,msh.c.Vc);
            e.t.a = B_2_1_1D.Update_et_a_f(s,pde.av);
            e.t.a = B_2_1_1D.Update_et_c  (msh,e.t.a);
            e.t.p = B_2_1_1D.Update_et_p_f(inp,msh,pde,s,stl,x.c);
            
            %  > Compute truncated terms (if requested).
            if ~inp.pa.adapt && inp.pl.tt
                s.nt  = inp.pl.nt;
                dfn_a = B_2_1_1D.Compute_dfA(s,stl,msh.f.Xv,pde.fn.f{1});
                e.t.a = A_2_1D.Compute_TTM  (inp,msh,s,stl,dfn_a);
            end
            %  > Update structures.
            [msh,pde] = B_2_1_1D.Set_struct(msh,pde,s,stl,x,e);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Initialize matrix coefficients.
        function [u,A,B] = Initialize_coeffs(j,NC,stl_s,sc,xf,A,B,bnd_s,bnd_v)
            %  > Initialize flags.
            u.c.uA = 0;
            u.w.uB = 0;
            u.e.uB = 0;
            
            if any(ismembc([j,j+1],stl_s))
                % >> A.
                %  > Update...
                A      = zeros(1,NC);
                u.c.uA = 1;
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
                        u.w.uB = 1;
                        u.w.fv = bnd_v(bnd_s == j);
                        u.w.Bf = xf{j}(l.w+1);
                    end
                    if ismembc(j+1,bnd_s)
                        u.e.uB = 1;
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
        function [s] = Update_s(inp,msh,pde,s,stl)
            %  > Auxiliary variables.
            NC = msh.c.NC;
            k  = 1:size(stl.s,2);
            Ac = zeros(NC);
            Bc = zeros(NC,1);
            
            %  > Update stencil.
            s = A_2_1D.Set_s(inp,msh,pde,s,stl);
            for i = 1:size(s.f,1)
                bnd_s{i} = find(~cellfun(@isempty,s.f(i,:)));
                bnd_v{i} = [s.bnd_v{i,bnd_s{i}}];
            end
            %  > Update A/B.
            for i = k
                if ~isempty(stl.s{i})
                    for j = 1:NC
                        %  > Initialize.
                        [u,A{i}(j,:),B{i}(j,1)] = ...
                            B_2_1_1D.Initialize_coeffs(j,NC,stl.s{i},s.c(i,:),s.xf(i,:),s.A{i}(j,:),s.B{i}(j,1),bnd_s{i},bnd_v{i});
                        
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
                else
                    A{i} = s.A{i};
                    B{i} = s.B{i};
                end
                %  > Add convective/diffusive contributions (cumulative matrices).
                Ac = Ac+s.vg(i).*A{i};
                Bc = Bc+s.vg(i).*B{i};
            end
            %  > Add source term.
            Bc   = Bc+pde.fn.st;
            %  > Update fields.
            s.A  = A;
            s.B  = B;
            s.Ac = Ac;
            s.Bc = Bc;
        end
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > 3.1.1. -------------------------------------------------------
        %  > Update 'pde.x.c' field (nodal solution).
        function [xc] = Update_xc(s)
            xc = s.Ac\s.Bc;
        end
        %  > 3.1.2. -------------------------------------------------------
        %  > Update 'pde.x.v' field (nodal values used to reconstruct face).
        function [vf] = Update_xv(s,x_c)
            %  > Auxiliary variables.
            [m,n] = size(s.xf);
            
            for i = 1:m
                for j = 1:n
                    if ~isempty(s.c{i,j})
                        %  > Nodal indices/values.
                        k = s.c{i,j};
                        v = x_c(k,1);
                        %  > Add boundary contribution?...
                        if ~isempty(s.f{i,j})
                            v = [v;s.bnd_v{i,j}];
                        end
                    else
                        v = s.bnd_v{i,j};
                    end
                    vf{i,j} = v;
                end
            end
        end
        %  > 3.1.3. -------------------------------------------------------
        %  > Update 'pde.x.f' field (reconstructed face values).
        function [xf] = Update_xf(s,x_c)
            %  > Auxiliary variables.
            [m,n] = size(s.xf);
            
            vf = B_2_1_1D.Update_xv(s,x_c);
            for i = 1:m
                for j = 1:n
                    xf(j,i) = s.xf{i,j}*vf{i,j};
                end
            end
        end
        % >> 3.2. ---------------------------------------------------------
        %  > 3.2.1. -------------------------------------------------------
        %  > Update 'pde.e.a.f' field (analytic face error).
        function [ef] = Update_ef(s,x_c,a_f)
            %  > Reconstructed face values (w/ input nodal field).
            x_ff          = B_2_1_1D.Update_xf(s,x_c);
            i             = 1:size(x_ff,2);
            %  > Error/absolute error distribution.
            ef.f    (:,i) = a_f(:,i)-x_ff(:,i);
            ef.f_abs(:,i) = abs(ef.f(:,i));
            %  > Mean error/absolute error.
            ef.n    (1,i) = mean(ef.f(:,i));
            ef.n_abs(1,i) = mean(ef.f_abs(:,i));
        end
        %  > 3.2.2. -------------------------------------------------------
        %  > Update 'pde.e.a.c' field (face cell error).
        function [ec] = Update_ec(x_c,a_c,Vc)
            %  > Error/absolute error distribution.
            ec.c    (:,1) = a_c(:,1)-x_c(:,1);
            ec.c_abs(:,1) = abs(ec.c);
            %  > Error/absolute error norms.
            ec_1    (:,1) = ec.c.*Vc;
            ec_1_abs(:,1) = ec.c_abs.*Vc;
            ec_2    (:,1) = ec_1.^2;
            ec_2_abs(:,1) = ec_1_abs.^2;
            ec.n    (1,1) = sum(ec_1)./sum(Vc);
            ec.n_abs(1,1) = sum(ec_1_abs)./sum(Vc);
            ec.n    (2,1) = sum(sqrt(ec_2))./sum(sqrt(Vc.^2));
            ec.n_abs(2,1) = sum(sqrt(ec_2_abs))./sum(sqrt(Vc.^2));
            ec.n    (3,1) = max(ec.c);
            ec.n_abs(3,1) = max(ec.c_abs);
        end
        %  > 3.2.3. -------------------------------------------------------
        %  > Update 'pde.e.t.a.f' field (analytic face truncation error).
        function [et] = Update_et_a_f(s,x)
            %  > Auxiliary variables.
            nc = length(s.vg);
            i  = 1:nc;
            j  = 1:nc+1;
            k  = j(end);

            %  > Reconstructed face values (w/ input nodal field).
            x.s             = B_2_1_1D.Update_xf(s,x.c);
            %  > (Weighted) truncation/absolute truncation error distribution.
            et.f      (:,i) = s.vg(i).*(x.f(:,i)-x.s(:,i));
            et.f      (:,k) = sum(et.f(:,i),2);
            et.f_abs  (:,j) = abs(et.f);
            %  > Mean/absolute mean truncation error.
            et.n.f    (:,j) = mean(et.f(:,j));
            et.n_abs.f(:,j) = mean(et.f_abs(:,j));
        end
        %  > 3.2.4. -------------------------------------------------------
        %  > Update 'pde.e.t.(...).c' field (cell truncation error).
        function [et] = Update_et_c(msh,et)
            %  > Truncation/absolute truncation error distribution.
            %  > Equivalent formulations (need to import x.c).
            %  et_c   (:,1) = s.Ac*av.c(:,1)-s.Bc(:,1);
            l               = 1:msh.c.NC;
            m               = l+1;
            k               = size(et.f,2);
            et.c      (l,1) = et.f(l,k)-et.f(m,k);
            et.c_abs  (:,1) = abs(et.c);
            
            %  > Error/absolute error norms.
            Vc              = msh.c.Vc;
            ec_1            = et.c.*Vc;
            ec_1_abs        = et.c_abs.*Vc;
            ec_2            = ec_1.^2;
            ec_2_abs        = ec_1_abs.^2;
            et.n.c    (1,:) = sum(ec_1,1)./sum(Vc);
            et.n.c    (2,:) = sum(sqrt(ec_2),1)./sum(sqrt(Vc.^2));
            et.n.c    (3,:) = max(et.c);
            et.n_abs.c(1,:) = sum(ec_1_abs,1)./sum(Vc);
            et.n_abs.c(2,:) = sum(sqrt(ec_2_abs),1)./sum(sqrt(Vc.^2));
            et.n_abs.c(3,:) = max(et.c_abs);
        end
        %  > 3.2.5. -------------------------------------------------------
        %  > Update 'pde.e.t.p' field (predicted/estimated error).
        function [et] = Update_et_p_f(inp,msh,pde,s,stl,x_c)
            %  > Auxiliary variables.
            ns = 2;
            
            for i = 1:ns+1
                if i == 1
                    sn{i} = s;
                else
                    for j = 1:size(stl.s,2)
                        k             = 2;
                        l             = j*size(stl.s,2)-1;
                        stl.p   (:,l) = inc(stl.p(:,l),k).i;
                        stl.s     {j} = transpose(1:msh.f.NF);
                    end
                    sn{i} = A_2_1D.Set_s(inp,msh,pde,s,stl);
                end
                fv_x{i} = B_2_1_1D.Update_xf(sn{i},x_c);
                if i ~= 1
                    for j = 1:size(stl.s,2)
                        et{i-1}.f(:,j) = sn{i}.vg(j).*(fv_x{i}(:,j)-fv_x{i-1}(:,j));
                    end
                    n                  = j+1;
                    et{i-1}.f    (:,n) = sum (et{i-1}.f,2);
                    et{i-1}.f_abs      = abs (et{i-1}.f);
                    et{i-1}.n.f        = mean(et{i-1}.f);
                    et{i-1}.n_abs.f    = mean(et{i-1}.f_abs);
                    et{i-1}            = B_2_1_1D.Update_et_c(msh,et{i-1});
                end
            end
        end
        
        %% > 4. -----------------------------------------------------------
        % >> 4.1. ---------------------------------------------------------
        %  > Compute derivatives (w/ analytic solution).
        function [df] = Compute_dfA(s,stl,Xv,f)
            %  > Auxiliary variables.
            l = size(stl.s,2);
            
            syms x;
            for i = 1:l
                for j = 0:s.nt(i)-1
                    k      = i*l-1;
                    n      = stl.p(k)+j+1;
                    dfn{i} = diff(f,x,n);
                    dfn{i} = matlabFunction(dfn{i});
                    if nargin(dfn{i}) ~= 0
                        df{i}(:,j+1) = dfn{i}(Xv)./factorial(n);
                    else
                        df{i}(:,j+1) = zeros(size(Xv));
                    end
                end
            end
        end
        % >> 4.2. ---------------------------------------------------------
        %  > Compute derivatives (w/ PDE solution).
        function [dfn] = Compute_dfN(s,tt,v)
            %  > Auxiliary variables.
            [m,n] = size(s.Inv);
            
            for i = 1:m
                for j = 1:n
                    k = 0;
                    for l = tt{i}
                        k           = k+1;
                        df          = zeros(1,size(s.Inv{i,j},1));
                        df    (1,l) = 1;
                        dfn{i}(j,k) = df*s.Inv{i,j}*v{i,j};
                    end
                end
            end
        end
        
        %% > 5. -----------------------------------------------------------
        %  > Update/set structure fields.
        function [msh,pde] = Set_struct(msh,pde,s,stl,x,e)
            %  > 'msh'.
            s.stl = stl;
            msh.s = s;
            msh   = Tools_1D.Order_msh(msh);
            %  > 'pde'.
            pde.x = x;
            pde.e = e;
            pde.e = Tools_1D.Order_pde_e(pde.e);
        end
    end
end
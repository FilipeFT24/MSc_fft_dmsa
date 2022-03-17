classdef B_2_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [msh,pde] = SetUp_p(inp,msh,pde,s,stl)
            %  > Update 's' structure and 'x' field.
            s     = B_2_1_1D.Update_s   (inp,msh,pde,s,stl);
            x.c   = B_2_1_1D.Update_xc  (s);
            %  > Compute cell/face/truncation error distribution/norms.
            e.a.c = B_2_1_1D.Update_ec_a(x.c,pde.av.c,msh.c.Vc);
            e.a.f = B_2_1_1D.Update_ef_a(x.c,pde.av.f,s);            
            e.a.t = B_2_1_1D.Update_et_a(pde.av,s,msh.c.Vc);
            e.p.t = B_2_1_1D.Update_et_p(inp,msh,pde,x.c,s,stl,2);
            e.p.c = B_2_1_1D.Update_ec_p(e.p.t,msh.c.Vc);

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
            Ac = zeros(NC);
            Bc = zeros(NC,1);
            
            %  > Update stencil.
            s = A_2_1D.Set_s(inp,msh,pde,s,stl);
            for i = 1:size(s.f,1)
                bnd_s{i} = find(~cellfun(@isempty,s.f(i,:)));
                bnd_v{i} = [s.bnd_v{i,bnd_s{i}}];
            end
            %  > Update A/B.
            for i = 1:size(stl.s,2)
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
                        k = s.c{i,j};
                        v = x_c(k,1);
                        %  > Add boundary contribution(?).
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
        %  > Update 'pde.e.a.c' field (cell error distribution/norms).
        function [ec] = Update_ec_a(xc,ac,Vc)
            %  > Error distribution.
            ec.c    (:,1) = ac(:,1)-xc(:,1);
            ec.c_abs(:,1) = abs(ec.c);
            %  > Error norms.
            ec.n          = LX.n(ec.c,Vc);
            ec.n_abs      = LX.n(ec.c_abs,Vc);
        end
        %  > 3.2.2. -------------------------------------------------------
        %  > Update 'pde.e.a.f' field (face error distribution/norms).
        function [ef] = Update_ef_a(xc,af,s)
            %  > Error distribution.
            i             = 1:size(af,2);
            xf            = B_2_1_1D.Update_xf(s,xc);
            ef.f    (:,i) = af(:,i)-xf(:,i);
            ef.f_abs(:,i) = abs(ef.f(:,i));
            %  > Error norms.
            ef.n    (:,i) = LX.n(ef.f(:,i));
            ef.n_abs(:,i) = LX.n(ef.f_abs(:,i));
        end
        %  > 3.2.3. -------------------------------------------------------
        %  > Update 'pde.e.a.t' field (truncation error distribution/norms).
        function [et] = Update_et_a(x,s,Vc)
            %  > Error distribution.
            [m,n]           = size(s.xf);
            x.s             = B_2_1_1D.Update_xf(s,x.c);
            i               = 1:m;
            j               = m+1;
            et.f      (:,i) = s.vg(i).*(x.f(:,i)-x.s(:,i));
            et.f      (:,j) = sum(et.f(:,i),2);
            et.f_abs        = abs(et.f);
            a               = 1:n-1;
            b               = a+1;
            et.c      (a,1) = et.f(a,j)-et.f(b,j);
            et.c_abs  (:,1) = abs(et.c);
            %  Remark: Equivalent formulation.
            %  et_c   (:,1)   = s.Ac*av.c(:,1)-s.Bc(:,1);            
            %  > Error norms.
            et.n.f          = LX.n(et.f);
            et.n_abs.f      = LX.n(et.f_abs);
            et.n.c          = LX.n(et.c,Vc);
            et.n_abs.c      = LX.n(et.c_abs,Vc);
        end
        %  > 3.2.4. -------------------------------------------------------
        %  > Update 'pde.e.p.t' field (predicted/estimated truncation error distribution/norms).
        function [et] = Update_et_p(inp,msh,pde,xc,s,stl,o)
            [m,n] = size(s.xf);
            %  > Stencil coefficients.
            for i = 1:o+1
                if i == 1
                    sn{i} = s;
                else
                    for j = 1:m
                        k          = 2;
                        l          = j*m-1;
                        stl.p(:,l) = stl.p(:,l)+k;
                        stl.s  {j} = transpose(1:n);
                    end
                    sn{i} = A_2_1D.Set_s(inp,msh,pde,s,stl);  
                end
                fv{i} = B_2_1_1D.Update_xf(sn{i},xc);
            end
            for i = 1:o
                %  > Error distribution.
                for j = 1:m
                    et{i}.f(:,j) = sn{i+1}.vg(j).*(fv{i+1}(:,j)-fv{i}(:,j));
                end
                k                = m+1;
                et{i}.f    (:,k) = sum(et{i}.f(:,1:m),2);
                et{i}.f_abs      = abs(et{i}.f);
                a                = 1:n-1; 
                b                = a+1;
                et{i}.c          = et{i}.f(a,k)-et{i}.f(b,k);
                et{i}.c_abs      = abs(et{i}.c);
                %  > Error norms.
                et{i}.n.f        = LX.n(et{i}.f);
                et{i}.n_abs.f    = LX.n(et{i}.f_abs);
                et{i}.n.c        = LX.n(et{i}.c,msh.c.Vc);
                et{i}.n_abs.c    = LX.n(et{i}.c_abs,msh.c.Vc);
                %  > Stencil coefficients.
                et{i}.s          = sn{i+1};
            end
        end
        %  > 3.2.5. -------------------------------------------------------
        %  > Update 'pde.e.p.c' field (predicted/estimated truncation error distribution/norms).
        function [ec] = Update_ec_p(et_p,Vc)
            for i = 1:size(et_p,2)
                %  > Error distribution.
                ec{i}.c       = inv(et_p{i}.s.Ac)*et_p{i}.c;
                ec{i}.c_abs   = abs(ec{i}.c);
                %  > Error norms.
                ec{i}.n.c     = LX.n(ec{i}.c,Vc);
                ec{i}.n_abs.c = LX.n(ec{i}.c_abs,Vc);
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
classdef B_2_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Update stencil coordinates,etc.
        function [s] = Update_1(inp,msh,f,s,u,add_b)
            s        = B_2_1_1D.Update_sc(inp,msh,f,s,u,add_b);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update stencil coefficients,etc.
        function [x] = Update_2(msh,f,s,u,x)
            x        = B_2_1_1D.Update_sx(msh,f,s,u,x);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Update matrices.
        function [m] = Update_3(msh,f,m,s,u,x)
            m        = B_2_1_1D.Update_mA(msh,f,m,s,u,x);
            m        = B_2_1_1D.Update_mB(msh,f,m,s,u,x);
        end
        % >> 1.4. ---------------------------------------------------------
        %  > Update nodal/face values, etc.
        function [x] = Update_4(msh,f,s,u,x)
            x        = B_2_1_1D.Update_xv   (msh,f,s,u,x);
            x        = B_2_1_1D.Update_cf   (u,x);
            x        = B_2_1_1D.Update_xf_xt(msh,f,s,u,x);
        end
        % >> 1.5. ---------------------------------------------------------
        %  > Update all.
        function [m,s,x] = Update_all(inp,msh,f,m,s,u,x,add_b,f_sol)
            s            = B_2_1_1D.Update_1 (inp,msh,f,s,u,add_b);
            x            = B_2_1_1D.Update_2 (msh,f,s,u,x);
            m            = B_2_1_1D.Update_3 (msh,f,m,s,u,x);
            if f_sol
                x.nv.x.c = B_2_1_1D.Update_xc(m.At,m.Bt);
            end
            x            = B_2_1_1D.Update_4 (msh,f,s,u,x);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Construct/update stencil coordinates.
        function [s] = Update_sc(inp,msh,f,s,u,add)
            %  > Auxiliary variables.
            Nf = msh.f.Nf;
            Xc = msh.c.Xc;
            Xv = msh.f.Xv;
            
            for i = 1:size(u.s,2)
                if isempty(u.s{i})
                    continue;
                else
                    for k  = u.s{i}'
                        %  > Stencil cell/face indices.
                        nb = ceil(u.p(k,i)./2);
                        ns = k-nb:k+nb-1;
                        %  > Check whether stencil reaches any of the boundaries and shift it accordingly.
                        if any(ns <= 0)
                            % >> WB.
                            nw = mcount(ns,0,'<');
                            ns = ns+nw;
                            sc = ns(ns ~= 0);
                            %  > Add 1 cell to the right(?).
                            if add
                                sc = [sc,sc(end)+1];
                            end
                            xt = [Xc(sc);Xv(1)]';
                        elseif any(ns >= Nf)
                            % >> EB.
                            ne = mcount(ns,Nf,'>');
                            ns = ns-ne;
                            sc = ns(ns ~= Nf);
                            %  > Add 1 cell to the left(?).
                            if add
                                sc = [sc(1)-1,sc];
                            end
                            xt = [Xc(sc);Xv(Nf)]';
                        else
                            sc = ns;
                            xt = Xc(sc)';
                        end
                        %  > Update 's' field.
                        s.c{k,i} = sc;
                        s.t{k,i} = sort(xt);
                    end
                end
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Construct/update stencil coefficients.
        function [x] = Update_sx(msh,f,s,u,x)
            for i = 1:size(u.s,2)
                if ~isempty(u.s{i})
                    for j = u.s{i}'
                        %  > Df.
                        xt        = s.t{j,i};
                        Df        = B_2_1_1D.Assemble_Df(f.bd,xt,msh.f.Xv(j),s.gv);
                        %  > Tf.
                        df        = zeros(1,length(xt));
                        df  (1,i) = 1;
                        Inv       = inv(Df);
                        Tf        = df*Inv;
                        %  > Update 'x' field.
                        x.Df{j,i} = Df;
                        x.if{j,i} = Inv;
                        x.Tf{j,i} = Tf;
                    end
                end
            end
        end
        %  > 2.2.1. -------------------------------------------------------
        %  > Check boundary type and assemble matrix Df accordingly.
        function [Df] = Assemble_Df(f_bd,xt,x,gv)
            %  > Auxiliary variables.
            lf = length(xt);
            Df = zeros (lf,lf);
            
            bd = ismembc(xt,f_bd.x);
            if any(bd) && f_bd.t(f_bd.x == xt(bd)) ~= "Dirichlet"
                switch f_bd.t(f_bd.x == xt(bd))
                    case "Neumann"
                        i        = 1:lf;
                        j        = 1:lf-1;
                        k        = xt == f_bd.x(f_bd.x == xt(bd));
                        Df(~k,i) = (xt(~k)-x)'.^(i-1);
                        Df( k,i) = [0,j.*(xt(k)-x)'.^(j-1)];
                    case "Robin"
                        i        = 1:lf;
                        j        = 1:lf-1;
                        k        = xt == f_bd.x(f_bd.x == xt(bd));
                        Df( i,i) = (xt(i)-x)'.^(i-1);
                        Df( k,i) = Df( k,i)+[0,j.*(xt(k)-x)'.^(j-1)].*gv;
                    otherwise
                        return;
                end
            else
                i       = 1:lf;
                Df(i,i) = (xt(i)-x)'.^(i-1);
            end
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > Construct/update matrix A.
        function [m] = Update_mA(msh,f,m,s,u,x)
            %  > Auxiliary variables.
            Nc = msh.c.Nc;
            
            % >> Set face equation(s) and construct/update A (cell dependent coefficient matrix).
            for i = 1:size(u.s,2)
                %  > Af (face contributions).
                if ~isempty(u.s{i})
                    for j = u.s{i}'
                        k                   = ismembc(s.t{j,i},f.bd.x);
                        m.Af{i}(j,:)        = zeros(1,Nc);
                        m.Af{i}(j,s.c{j,i}) = x.Tf{j,i}(~k);
                    end
                    %  > Cell indices to be updated.
                    a{i}(:,1) = unique([msh.f.c{u.s{i}}]);
                    
                    %  > Ac (cell contributions).
                    if ~isempty(a{i})
                        for j = a{i}'
                            %  > Add west(w)/east(e) face contributions...
                            k            = msh.c.f.f{j};
                            m.Ac{i}(j,:) = zeros(1,Nc);
                            for n = 1:length(k)
                                m.Ac{i}(j,:) = m.Ac{i}(j,:)+m.Af{i}(k(n),:).*msh.c.f.Nf{j}(n);
                            end
                        end
                    end
                end
                m.nnz.Ac(i) = B_2_1_1D.Set_nnz(m.Ac{i});
                m.nnz.Af(i) = B_2_1_1D.Set_nnz(m.Af{i});
            end
            %  > Cell(s) to be updated.
            b = unique(cat(1,a{:}));
            
            %  > At (cumulative/total matrix).
            if ~isempty(b)
                for i = b'
                    %  > Add convective/diffusive contributions...
                    m.At(i,:) = zeros(1,Nc);
                    for j = 1:size(u.s,2)
                        m.At(i,:) = m.At(i,:)+s.v(j).*m.Ac{j}(i,:);
                    end
                end
            end
            m.nnz.At = B_2_1_1D.Set_nnz(m.At);
        end
        % >> 3.2. ---------------------------------------------------------
        %  > Construct/update matrix B.
        function [m] = Update_mB(msh,f,m,s,u,x)
            %  > Auxiliary variables.
            Nf   = msh.f.Nf;
            k_bf = zeros(Nf,size(u.s,2));
            
            % >> Set face equation(s) and construct/update B (face dependent coefficient matrix w/ source term contribution).
            for i = 1:size(u.s,2)
                %  > Bf (face contributions).
                if ~isempty(u.s{i})
                    for j = u.s{i}'
                        k = ismembc(s.t{j,i},f.bd.x);
                        if any(k)
                            k_bf   (j,i) = 1;
                            m.Bf{i}(j,1) = x.Tf{j,i}(k).*f.bd.v(f.bd.x == s.t{j,i}(k));
                        end
                    end
                end
                %  > Cell indices to be updated.
                a{i}(:,1) = find(k_bf(:,i));
                b{i}(:,1) = unique([msh.f.c{a{i}(:,1)}]);
                
                %  > Bc (cell contributions).
                if ~isempty(b{i})
                    for j = b{i}'
                        %  > Add west(w)/east(e) face contributions...
                        k            = msh.c.f.f{j};
                        m.Bc{i}(j,1) = 0;
                        for l = 1:length(k)
                            m.Bc{i}(j,1) = m.Bc{i}(j,1)-m.Bf{i}(k(l),1).*msh.c.f.Nf{j}(l);
                        end
                    end
                end
            end
            %  > Cell(s) to be updated.
            c = unique(cat(1,b{:}));
            
            %  > Bt (cumulative/total matrix).
            if ~isempty(c)
                for i = c'
                    %  > Add convective/diffusive contributions...
                    m.Bt(i,1) = f.fv(i);
                    for j = 1:size(u.s,2)
                        m.Bt(i,1) = m.Bt(i,1)+s.v(j).*m.Bc{j}(i,1);
                    end
                end
            end
        end
        % >> 3.3. ---------------------------------------------------------
        %  > Compute nnz of each matrix.
        function [nnz_m] = Set_nnz(m)
            nnz_m = nnz(m);
        end
        
        %% > 4. -----------------------------------------------------------
        % >> 4.1. ---------------------------------------------------------
        %  > Update 'x.nv.x.c' field (nodal solution).
        function [xc] = Update_xc(At,Bt)
            xc = At\(Bt);
        end
        % >> 4.2. ---------------------------------------------------------
        %  > Update 'x.vf' field (nodal values used to fit face polynomial).
        function [x] = Update_xv(msh,f,s,u,x)
            for i = ["a","x"]
                for j = 1:size(u.s,2)
                    if ~isempty(u.s{j})
                        for k = u.s{j}'
                            l = ismembc(s.t{k,j},f.bd.x);
                            %  > Boundary face contribution(s).
                            if any(l)
                                x.vf.(i){k,j}(l,1) = f.bd.v(f.bd.x == s.t{k,j}(l));
                            end
                            %  > Cell contribution(s).
                            x.vf.(i){k,j}(~l,1) = x.nv.(i).c(s.c{k,j});
                        end
                    end
                end
            end
        end
        % >> 4.3. ---------------------------------------------------------
        %  > Update 'x.cf' field (fitted polynomial coefficients).
        function [x] = Update_cf(u,x)
            for i = ["a","x"]
                for j = 1:size(u.s,2)
                    if ~isempty(u.s{j})
                        for k = u.s{j}'
                            x.cf.(i){k,j} = x.if{k,j}*x.vf.(i){k,j};
                        end
                    end
                end
            end
        end
        % >> 4.4. ---------------------------------------------------------
        %  > Update 'x.f' and 'x.t' fields (nodal face/total values).
        function [x] = Update_xf_xt(msh,f,s,u,x)
            for i = ["a","x"]
                for j = 1:size(u.s,2)
                    if ~isempty(u.s{j})
                        for k = u.s{j}'
                            %  > 'xf'.
                            x.xf.(i)(k,j) = x.Tf{k,j}*x.vf.(i){k,j};
                            %  > 'xt'.
                            n_cf          = length(x.cf.(i){k,j});
                            D_cf          = zeros(n_cf);
                            x.xt.(i){k,j} = zeros(n_cf,1);
                            x_ff          = msh.f.Xv(k);
                            t             = s.t{k,j};
                            %  > Column #1: Interpolated value.
                            %  > Column #2: Absolute error of interpolated value.
                            %           └─ x.xt.(i){k,j}(:,2) = abs(f.fh{j}(t)'-x.xt.(i){k,j}(:,1));
                            switch j
                                case 1
                                    l                  = 1:n_cf;
                                    x.xt.(i){k,j}(:,1) = (x_ff-t').^(l-1)*x.cf.(i){k,j}(l);
                                case 2
                                    l                  = 1:n_cf-1;
                                    m                  = l+1;
                                    D_cf         (:,m) = l.*(x_ff-t').^(l-1);
                                    x.xt.(i){k,j}(:,1) = D_cf*x.cf.(i){k,j};
                            end
                        end
                    end
                end
            end
        end
        
        %% > 5. -----------------------------------------------------------
        % >> 5.1. ---------------------------------------------------------
        %  > Update 'e.p(...)' field (predicted/estimated cell/face truncation error distribution/norms).
        function [ep] = Update_ep(ep,mo,s,xo,xn,Vc)
            %  > \tau_f(\phi) & \tau(\nabla\phi).
            for i = 1:size(ep.t.f,2)-1
                ep.t.f(:,i) = s.v(i).*(xn.xf.x(:,i)-xo.xf.x(:,i));
            end
            %  > Update remaining error fields...
            ep = Tools_1D.Set_e(ep,mo,Vc);
        end
        % >> 5.2. ---------------------------------------------------------
        %  > Update 'e.a(...)' field (analytic cell/face truncation error distribution/norms).
        function [ea] = Update_ea(ea,m,s,x,Vc)
            %  > \tau_f(\phi) & \tau(\nabla\phi).
            for i = 1:size(ea.t.f,2)-1
                ea.t.f(:,i) = s.v(i).*(x.nv.a.f(:,i)-x.xf.a(:,i));
            end
            %  > Update remaining error fields...
            ea = Tools_1D.Set_e(ea,m,Vc);
        end
        % >> 5.3. ---------------------------------------------------------
        %  > Update 'e.d(...)' field (difference analytic/predicted).
        function [ed] = Update_ed(ea,ed,ep,Vc)
            %  > Field.
            fld       = ["c","t"];
            g{1}(1)   = ["c"];
            g{2}(1,:) = ["c","f"];
            %  > Error distribution.
            for i = 1:size(g,2)
                a = fld(i);
                for j = 1:size(g{i},2)
                    b          = g{i}(j);
                    c          = strcat(b,"_abs");
                    ed.(a).(b) = ea.(a).(b)-ep.(a).(b);
                    ed.(a).(c) = abs(ed.(a).(b));
                end
            end
            %  > Error norms.
            ed.c.n       = Tools_1D.Set_n(ed.c.c,Vc);
            ed.c.n_abs   = Tools_1D.Set_n(ed.c.c_abs,Vc);
            ed.t.n.c     = Tools_1D.Set_n(ed.t.c,Vc);
            ed.t.n.f     = Tools_1D.Set_n(ed.t.f);
            ed.t.n_abs.c = Tools_1D.Set_n(ed.t.c_abs,Vc);
            ed.t.n_abs.f = Tools_1D.Set_n(ed.t.f_abs);
        end
        % >> 5.4 ---------------------------------------------------------
        %  > Update 'e.(...)' field (error).
        function [e,x_s] = Update_e(inp,msh,e,f,m,s,u,x,add_b)
            %  > Auxiliary variables.
            f_sol  = 0;
            Nf     = msh.f.Nf;
            Vc     = msh.c.Vc;
            x_s{1} = x;
            
            for i = 1:inp.pa.ns
                %  > Update field 'u'.
                u       = B_2_1_1D.Set_upd_p (u,Nf);
                %  > Assign fields 'e.a', 'mp' and 'xp'.
                e.a{i}  = B_2_1_1D.Update_ea (e.a{i},m,s,x,Vc);
                mp      = m;
                xp      = x;
                %  > Update stencil...
                [m,s,x] = B_2_1_1D.Update_all(inp,msh,f,m,s,u,x,add_b,f_sol);
                %  > Assign fields 'x', 'e.d' and 'e.p'.
                j       = i+1;
                x_s{j}  = x;
                e.p{i}  = B_2_1_1D.Update_ep (e.p{i},mp,s,xp,x,Vc);
                e.d{i}  = B_2_1_1D.Update_ed (e.a{i},e.d{i},e.p{i},Vc);
            end
        end
        %  > 5.4.1. -------------------------------------------------------
        %  > Auxiliary function (increase method's order).
        function [u] = Set_upd_p(u,Nf)
            A = 2;
            for i = 1:size(u.s,2)
                if ~isempty(u.s{i})
                    for j = u.s{i}'
                        u.p(j,i) = u.p(j,i)+A;
                    end
                end
                u_s{i}(:,1) = 1:Nf;
                u.s{i}      = u_s{i};
            end
        end
    end
end
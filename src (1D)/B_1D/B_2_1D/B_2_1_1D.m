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
        function [x] = Update_4(f,s,u,x)           
            x        = B_2_1_1D.Update_xv(f,s,u,x);  
            x        = B_2_1_1D.Update_cf(u,x);
            x        = B_2_1_1D.Update_xf(u,x);
        end
        % >> 1.5. ---------------------------------------------------------
        %  > Update all.
        function [m,s,x] = Update_all(inp,msh,f,m,s,u,x,add_b,upd)
            %  > Update stencil.
            s            = B_2_1_1D.Update_1 (inp,msh,f,s,u,add_b);
            x            = B_2_1_1D.Update_2 (msh,f,s,u,x);
            %  > Update matrices(?).
            if upd(1)
                m        = B_2_1_1D.Update_3 (msh,f,m,s,u,x);
            end
            %  > Update nodal solution(?).
            if upd(2)
                x.nv.x.c = B_2_1_1D.Update_xc(m.At,m.Bt);
            end
            %  > Update face values(?).
            if upd(3)
                x        = B_2_1_1D.Update_4 (f,s,u,x);
            end
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
                        p  = u.p(k,i);
                        %  > Stencil cell/face indices.
                        nb = ceil(p./2);
                        ns = k-nb:k+nb-1;
                        %  > Check whether stencil reaches any of the boundaries and shift it accordingly.
                        if any(ns <= 1)
                            % >> WB.
                            nw = mcount(ns,0,'<');
                            ns = ns+nw;
                            %  > Add cell(s) to the right(?).
                            if ~add
                                sc = ns(ns ~= 0);
                                xt = Xc(sc)';
                                if any(ns == 0)
                                    xt = [xt,f.bd.x(1)];
                                end
                            else
                                
                            
                                
                                
                                
                                sc = 1:p+1;
                                xt = [Xc(sc)',f.bd.x(1)];
                              
  
                                
                                
                                
                                
                                
                            end                         
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
        function [x] = Update_xv(f,s,u,x)
            for i = u.f
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
            for i = u.f
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
        %  > Update 'x.f' field (nodal face values).
        function [x] = Update_xf(u,x)
            for i = u.f
                for j = 1:size(u.s,2)
                    if ~isempty(u.s{j})
                        for k = u.s{j}'
                            x.xf.(i)(k,j) = x.Tf{k,j}*x.vf.(i){k,j};
                        end
                    end
                end
            end
        end
        
        %% > 5. -----------------------------------------------------------
        % >> 5.1. ---------------------------------------------------------
        %  > Update 'e.p(...)' field (predicted/estimated cell/face truncation error distribution/norms).
        function [ed,ep,x] = Update_ed_ep(eda,ed,ep,f,m,s,u,x,Vc,fs)
            %  > \tau_f(\phi), \tau(\nabla\phi), \tau_f and \tau_c.
            ep = Tools_1D.Set_1_e(ep,s{1}.v,x{1}.xf.x,x{2}.xf.x);
            %  > Update remaining error fields...
            if ~fs(2)
                if ~fs(1)
                    %  > w/ LO (j=2).
                    j = 2;
                else
                    %  > w/ HO (j=1).
                    j = 1;
                end
                ep = Tools_1D.Set_2_e(ep,m{j},Vc);
            else
                %  > Add as source(?).
                j  = 1;
                x  = B_2_1_1D.Add_ep_tc(ep,f,m{j},s,u,x);
                ep = Tools_1D.Set_1_e  (ep,s{1}.v,x{1}.xf.x,x{2}.xf.x);
                ep = Tools_1D.Set_2_e  (ep,m{j},Vc);
            end
            ed = B_2_1_1D.Update_ed(eda,ed,ep,Vc);
        end
        %  > 5.1.1. -------------------------------------------------------
        %  > Auxiliary function (add \tau_c(p) as source on the RHS).  
        function [x] = Add_ep_tc(ep,f,m,s,u,x)
            %  > Update cell values...
            xc = B_2_1_1D.Update_xc(m.At,m.Bt+ep.t.c);
            %  > Update "x.(...).x" field only.
            for i = 1:size(x,2)
                u{i}.f      = "x";   
                x{i}.nv.x.c = xc;
                x{i}        = B_2_1_1D.Update_4(f,s{i},u{i},x{i});
            end
        end
        % >> 5.2. ---------------------------------------------------------
        %  > Update 'e.a(...)' field (analytic cell/face truncation error distribution/norms).
        function [ea] = Update_ea(ea,m,v,x,Vc)
            %  > \tau_f(\phi), \tau(\nabla\phi), \tau_f and \tau_c.
            ea = Tools_1D.Set_1_e(ea,v,x.xf.a,x.nv.a.f);
            %  > Update remaining error fields...
            ea = Tools_1D.Set_2_e(ea,m,Vc);
        end
        % >> 5.3. ---------------------------------------------------------
        %  > Update 'e.da(...)' field (analytic cell/face (difference) truncation error distribution/norms).
        function [eda] = Update_eda(eda,m,v,x,Vc)
            %  > \tau_f(\phi), \tau(\nabla\phi), \tau_f and \tau_c.
            eda = Tools_1D.Set_1_e(eda,v,x{1}.xf.a,x{2}.xf.a);
            %  > Update remaining error fields...
            eda = Tools_1D.Set_2_e(eda,m,Vc);
        end
        % >> 5.4. ---------------------------------------------------------
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
        % >> 5.5 ---------------------------------------------------------
        %  > Update 'e.(...)' field (error).
        function [e,xs] = Update_e(inp,msh,e,f,m,s,u,x,add_b)
            %  > Auxiliary variables.
            Nf    = msh.f.Nf;
            Vc    = msh.c.Vc;
            v     = s.v;
            %  > Options.
            fs(1) = 1;   %  > Use higher-order solution(?).
            fs(2) = 0;   %  > Add lower-order (predicted) cell truncation error as source term(?).
            if all(fs)   %  > Allow only w/ lower-order solution.
                return;
            end
                        
            %  > Assign...
            if ~fs(1)
                x.nv.x.c = B_2_1_1D.Update_xc(m.At,m.Bt);
                x        = B_2_1_1D.Update_4 (f,s,u,x);
            end
            ms{1} = m;
            ss{1} = s;
            us{1} = u;
            xs{1} = x;
            ns    = inp.pa.ns;
            %  > #1: Update stencil...
            for i = 1:ns
                %  > Update fields 'm', 's', 'u' and 'x'.
                u       = B_2_1_1D.Set_upd_p (u,Nf);
                [m,s,x] = B_2_1_1D.Update_all(inp,msh,f,m,s,u,x,add_b,[1,fs(1),1]);
                %  > Assign fields 'ms', 'ss', 'us' and 'xs'.               
                if fs(1)
                    xs{i}.nv.x.c = x.nv.x.c;
                    xs{i}        = B_2_1_1D.Update_4(f,ss{i},us{i},xs{i});
                end
                ms{i+1} = m;
                ss{i+1} = s;
                us{i+1} = u;
                xs{i+1} = x;
            end
            %  > #2: Update fields 'e.a' and 'e.da': if the (predicted) cell truncation error is added as a source term, no need to update matrices, since e=A(LO)\(\tau_c).
            for i = 1:ns+1
                e.a{i} = B_2_1_1D.Update_ea(e.a{i},ms{i},v,xs{i},Vc);
                if i ~= 1
                    j         = i-1:i;
                    e.da{i-1} = B_2_1_1D.Update_eda(e.da{i-1},ms{i},v,xs(j),Vc);
                end
            end
            %  > #3: Update fields 'e.d', 'e.p' and 'xs'.
            for i = 1:ns
                j                     = i:i+1;
                [e.d{i},e.p{i},xs(j)] = ...
                    B_2_1_1D.Update_ed_ep(e.da{i},e.d{i},e.p{i},f,ms(j),ss(j),us(j),xs(j),Vc,fs);
            end
        end
        %  > 5.5.1. -------------------------------------------------------
        %  > Auxiliary function (increase method's order).
        function [u] = Set_upd_p(u,Nf)
            A = 2;
            for i = 1:size(u.s,2)
                if ~isempty(u.s{i})
                    u.p(u.s{i},i) = u.p(u.s{i},i)+A;       
                end
                u_s{i}(:,1) = 1:Nf;
                u.s{i}      = u_s{i};
            end
        end
    end
end
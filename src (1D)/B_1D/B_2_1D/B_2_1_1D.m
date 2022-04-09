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
        function [x] = Update_2(msh,s,u,x)
            x        = B_2_1_1D.Update_sx(msh,s,u,x);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Update matrices.
        function [m] = Update_3(msh,f,m,s,u,x)
            m        = B_2_1_1D.Update_mA(msh,m,s,u,x);
            m        = B_2_1_1D.Update_mB(msh,f,m,s,u,x);
        end
        % >> 1.4. ---------------------------------------------------------
        %  > Update nodal/face values, etc.
        function [x] = Update_4(f,s,u,x)
            x        = B_2_1_1D.Update_xv(s,u,x);
            x        = B_2_1_1D.Update_xf(f,s,u,x);
        end
        % >> 1.5. ---------------------------------------------------------
        %  > Update all.
        function [m,s,x] = Update_all(inp,msh,f,m,s,u,x,add_b,f_sol)
            s            = B_2_1_1D.Update_1 (inp,msh,f,s,u,add_b);
            x            = B_2_1_1D.Update_2 (msh,s,u,x);
            m            = B_2_1_1D.Update_3 (msh,f,m,s,u,x);
            if f_sol
                x.nv.x.c = B_2_1_1D.Update_xc(m.At,m.Bt);
            end
            x            = B_2_1_1D.Update_4 (f,s,u,x);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Construct/update stencil coordinates.
        function [s] = Update_sc(inp,msh,f,s,u,add)
            %  > Auxiliary variables.
            Xc = msh.c.Xc;
            Xv = msh.f.Xv;
            Nf = msh.f.Nf;
            nv = size(u.s,2);
            
            for i = 1:nv
                if isempty(u.s{i})
                    continue;
                else
                    j = i.*nv-1;
                    k = i.*nv;
                    for l  = 1:size(u.s{i},1)
                        %  > Face index.
                        m  = u.s{i}(l);
                        %  > Check polynomial order (odd/even).
                        nb = ceil(u.p(m,j)./2);
                        nl = -nb;
                        nr =  nb;
                        if ~rem(u.p(m,i),2)
                            if u.p(m,k) < 0
                                %  > UDS.
                                nl = nl+u.p(m,k);
                            else
                                %  > DDS.
                                nr = nr+u.p(m,k);
                            end
                        end
                        ns = m+nl:m+nr-1;
                        
                        %  > Check whether stencil reaches any of the boundaries and shift it accordingly.
                        if any(ns <= 0)
                            % >> WB.
                            %  > Stencil indices.
                            nw = mcount(ns,0,'<');
                            ns = ns+nw;
                            sc = ns(ns ~= 0);
                            %  > Add 1 cell to the right.
                            if add
                                sc = [sc,sc(end)+1];
                            end
                            %  > Boundary type/value.
                            bt = inp.pv.b(1);
                            bv = B_2_1_1D.Set_bnd_v(inp,f.av.f,1,bt);
                            %  > Stencil coordinates.
                            xt = [Xc(sc);Xv(1)]';
                        elseif any(ns >= Nf)
                            % >> WB.
                            %  > Stencil indices.
                            ne = mcount(ns,Nf,'>');
                            ns = ns-ne;
                            sc = ns(ns ~= Nf);
                            %  > Add 1 cell to the left.
                            if add
                                sc = [sc(1)-1,sc];
                            end
                            %  > Boundary type/value.
                            bt = inp.pv.b(2);
                            bv = B_2_1_1D.Set_bnd_v(inp,f.av.f,Nf,bt);
                            %  > Stencil coordinates.
                            xt = [Xc(sc);Xv(Nf)]';
                        else
                            %  > Stencil indices.
                            sc = ns;
                            bt = [];
                            bv = [];
                            %  > Stencil coordinates.
                            xt = Xc(sc)';
                        end
                        %  > Update 's' field.
                        s.c {m,i} = sc;
                        s.t {m,i} = xt;
                        s.bt{m,i} = bt;
                        s.bv{m,i} = bv;
                    end
                end
            end
        end
        %  > 2.1.1. -------------------------------------------------------
        %  > Set boundary values.
        function [b_v] = Set_bnd_v(inp,func,f,b_t)
            switch string(b_t)
                case "Dirichlet"
                    b_v = func(f,1);
                case "Neumann"
                    b_v = func(f,2);
                case "Robin"
                    g_v = inp.pv.v(2)./inp.pv.v(1);
                    b_v = func(f,1)+g_v.*func(f,2);
                otherwise
                    b_v = [];
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Construct/update stencil coefficients.
        function [x] = Update_sx(msh,s,u,x)
            %  > Auxiliary variables.
            Xv = msh.f.Xv;
            
            for i = 1:size(u.s,2)
                if ~isempty(u.s{i})
                    for j = 1:size(u.s{i},1)
                        %  > Face index.
                        a     = u.s{i}(j);
                        %  > Df.
                        b.f   = Xv  (a);
                        b.xt  = s.t {a,i};
                        b.bt  = s.bt{a,i};
                        b.g_v = abs(s.v(2)./s.v(1));
                        Df    = B_2_1_1D.Assemble_Df(b);
                        %  > cf.
                        if length(b.xt) == 1
                            if i == 2
                                warndlg('1st order UDS/DDS cannot be used to discretize the diffusive term.');
                                break;
                            else
                                %                                %  > remove...
                                %                                 df  = 1;
                                %                                 Inv = inv(Df);
                                %                                 Tf  = df*Inv;
                                %                                 cf  = Tf(n,:);
                            end
                        else
                            df      = zeros(1,length(b.xt));
                            df(1,i) = 1;
                            Inv     = inv(Df);
                            Tf      = df*Inv;
                        end
                        %  > Update 'x' field.
                        x.Tf{a,i}   = Tf;
                        x.if{a,i}   = Inv;
                    end
                end
            end
        end
        %  > 2.2.1. -------------------------------------------------------
        %  > Assemble matrix Df.
        function [Df] = Assemble_Df(b)
            %  > Auxiliary variables.
            len = length(b.xt);
            j   = 1:len;
            Df  = zeros(len);
            
            switch string(b.bt)
                case "Dirichlet"
                    Df(j,j) = (b.xt(j)-b.f)'.^(j-1);
                case "Neumann"
                    k       = 1:len-1;
                    Df(k,j) = (b.xt(k)-b.f)'.^(j-1);
                    l       = len;
                    m       = k+1;
                    Df(l,1) = 0;
                    Df(l,m) = k.*(b.xt(l)-b.f).^(k-1);
                case "Robin"
                    g_v     = b.g_v;
                    k       = 1:len-1;
                    Df(k,j) = (b.xt(k)-b.f)'.^(j-1);
                    l       = len;
                    lv  (j) = (b.xt(l)-b.f)'.^(j-1);
                    m       = k+1;
                    lg  (1) = 0;
                    lg  (m) = k.*(b.xt(l)-b.f).^(k-1);
                    Df(l,j) = lv(j)+g_v*lg(j);
                otherwise
                    Df(j,j) = (b.xt(j)-b.f)'.^(j-1);
            end
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > Construct/update matrix A.
        function [m] = Update_mA(msh,m,s,u,x)
            %  > Auxiliary variables.
            Nc = msh.c.Nc;
            
            % >> Set face equation(s) and construct/update A (cell dependent coefficient matrix).
            for i = 1:size(u.s,2)
                %  > Af (face contributions).
                if ~isempty(u.s{i})
                    for j = 1:size(u.s{i},1)
                        %  > Face index/cell indices used to reconstruct face.
                        a            = u.s{i}(j);
                        b            = s.c{a,i};
                        %  > Overwrite/update field...
                        m.Af{i}(a,:) = zeros(1,Nc);
                        m.Af{i}(a,b) = x.Tf{a,i}(1:length(b));
                    end
                    %  > Cell indices to be updated.
                    c{i}(:,1) = unique([msh.f.c{u.s{i}}]);
                    
                    %  > Ac (cell contributions).
                    if ~isempty(c{i})
                        for j = 1:size(c{i},1)
                            %  > Cell 'j' face indices.
                            k = c{i}(j);
                            l = msh.c.f.f{k};
                            %  > Overwrite/update field and add west(w)/east(e) face contributions...
                            m.Ac{i}(k,:) = zeros(1,Nc);
                            for n = 1:length(l)
                                m.Ac{i}(k,:) = m.Ac{i}(k,:)+m.Af{i}(l(n),:).*msh.c.f.Nf{k}(n);
                            end
                        end
                    end
                end
                m.nnz.Ac(i) = B_2_1_1D.Set_nnz(m.Ac{i});
                m.nnz.Af(i) = B_2_1_1D.Set_nnz(m.Af{i});
            end
            %  > Cell(s) to be updated.
            d = unique(cat(1,c{:}));
            
            %  > At (cumulative/total matrix).
            if ~isempty(d)
                for i = 1:size(d,1)
                    %  > Overwrite/update field...
                    e         = d(i);
                    m.At(e,:) = zeros(1,Nc);
                    %  > Add convective/diffusive contributions...
                    for j = 1:size(u.s,2)
                        m.At(e,:) = m.At(e,:)+s.v(j).*m.Ac{j}(e,:);
                    end
                end
            end
            m.nnz.At = B_2_1_1D.Set_nnz(m.At);
        end
        % >> 3.2. ---------------------------------------------------------
        %  > Construct/update matrix B.
        function [m] = Update_mB(msh,f,m,s,u,x)
            % >> Set face equation(s) and construct/update B (face dependent coefficient matrix w/ source term contribution).
            for i = 1:size(u.s,2)
                %  > Bf (face contributions).
                if ~isempty(u.s{i})
                    for j = 1:size(u.s{i},1)
                        %  > Face index.
                        a = u.s{i}(j);
                        if ~isempty(s.bv{a,i})
                            %  > Overwrite/update field...
                            m.Bf{i}(j,1) = x.Tf{a,i}(end).*s.bv{a,i};
                        end
                    end
                end
                %  > Cell indices to be updated.
                b{i}(:,1) = find(~cellfun(@isempty,s.bv(:,i)));
                c{i}(:,1) = unique([msh.f.c{b{i}(:,1)}]);
                
                %  > Bc (cell contributions).
                for j = 1:size(c{i},1)
                    if ~isempty(c{i})
                        %  > Cell 'j' face indices.
                        k = c{i}(j);
                        l = msh.c.f.f{k};
                        
                        %  > Overwrite/update field and add west(w)/east(e) face contributions...
                        m.Bc{i}(k,1) = 0;
                        for n = 1:length(l)
                            m.Bc{i}(k,1) = m.Bc{i}(k,1)-m.Bf{i}(l(n),1).*msh.c.f.Nf{k}(n);
                        end
                    end
                end
            end
            %  > Cell(s) to be updated.
            d = unique(cat(1,c{:}));
            
            %  > Bt (cumulative/total matrix).
            if ~isempty(d)
                for i = 1:size(d,1)
                    %  > Overwrite/update field...
                    e         = d(i);
                    m.Bt(e,1) = f.fv(e);
                    %  > Add convective/diffusive contributions...
                    for j = 1:size(u.s,2)
                        m.Bt(e,1) = m.Bt(e,1)+s.v(j).*m.Bc{j}(e,1);
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
        function [d] = Update_xc(At,Bt)
            d = At\(Bt);
        end
        % >> 4.2. ---------------------------------------------------------
        %  > Update 'x.vf' field (nodal values used to fit face polynomial).
        function [x] = Update_xv(s,u,x)
            for i = ["a","x"]
                for j = 1:size(u.s,2)
                    if ~isempty(u.s{j})
                        for k = 1:size(u.s{j},1)
                            l = u.s{j}(k);
                            %  > Overwrite/update field...
                            if ~isempty(s.c{l,j})
                                %  > Cell contribution(s).
                                m     = s.c{l,j};
                                v.(i) = x.nv.(i).c(m,1);
                                %  > Boundary contribution(s).
                                if ~isempty(s.bv{l,j})
                                    v.(i) = [v.(i);s.bv{l,j}];
                                end
                            else
                                v.(i) = s.bv{l,j};
                            end
                            x.vf.(i){l,j} = v.(i);
                        end
                    end
                end
            end
        end
        % >> 4.3. ---------------------------------------------------------
        %  > Update 'x.cf' and 'x.xf' (fitted face coefficients/values).
        function [x] = Update_xf(f,s,u,x)
            for i = ["a","x"]
                for j = 1:size(u.s,2)
                    if ~isempty(u.s{j})
                        for k = 1:size(u.s{j},1)
                            l             = u.s{j}(k);
                            x.xf.(i)(l,j) = x.Tf{l,j}*x.vf.(i){l,j};
                            x.cf.(i){l,j} = x.if{l,j}*x.vf.(i){l,j};
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
        %  > Update 'e.d(...)' field.
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
        %  > Update 'e.(...)' field.
        function [e] = Update_e(inp,msh,e,f,m,s,u,x,add_b)
            %  > Auxiliary variables.
            f_sol = 0;
            Nf    = msh.f.Nf;
            Vc    = msh.c.Vc;
            
            for i = 1:inp.pa.ns
                %  > Update field 'u'.
                u       = B_2_1_1D.Set_upd_p (u,Nf);
                %  > Assign fields 'e.a', 'mp' and 'xp'.
                e.a{i}  = B_2_1_1D.Update_ea (e.a{i},m,s,x,Vc);
                mp      = m;
                xp      = x;
                %  > Update stencil...
                [m,s,x] = B_2_1_1D.Update_all(inp,msh,f,m,s,u,x,add_b,f_sol);
                %  > Assign fields 'e.d' and 'e.p'.
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
                    for j = 1:size(u.s{i},1)
                        k        = u.s{i}(j);
                        l        = i*size(u.s,2)-1;
                        u.p(k,l) = u.p(k,l)+A;
                    end
                end
                u_s{i}(:,1) = 1:Nf;
                u.s{i}      = u_s{i};
            end
        end
    end
end
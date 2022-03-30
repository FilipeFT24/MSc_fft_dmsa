classdef B_2_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Update stencil coordinates,etc.
        function [s] = Update_1(inp,msh,pde,s,u,add)
            s        = B_2_1_1D.Update_sc(inp,msh,pde,s,u,add);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update stencil coefficients,etc.
        function [x] = Update_2(inp,msh,s,u,x)
            x        = B_2_1_1D.Update_sx(inp,msh,s,u,x);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Update matrices.
        function [m] = Update_3(msh,pde,m,s,u,x)
            m        = B_2_1_1D.Update_mA(msh,m,s,u,x);
            m        = B_2_1_1D.Update_mB(msh,pde,m,s,u,x);
        end
        % >> 1.4. ---------------------------------------------------------
        %  > Update nodal/face values.
        function [x] = Update_4(s,u,x)
            x        = B_2_1_1D.Update_xv(s,u,x);
            x        = B_2_1_1D.Update_xf(u,x);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Construct/update stencil coordinates.
        function [s] = Update_sc(inp,msh,pde,s,u,add)
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
                            bv = B_2_1_1D.Set_bnd_v(inp,pde.av.f,1,bt);
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
                            bv = B_2_1_1D.Set_bnd_v(inp,pde.av.f,Nf,bt);
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
        function [x] = Update_sx(inp,msh,s,u,x)
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
                        b.g_v = inp.pv.v(2)./inp.pv.v(1);
                        Df    = B_2_1_1D.Assemble_Df(b);
                        %  > cf.
                        if length(b.xt) == 1
                            if i == 2
                                warndlg('1st order UDS/DDS cannot be used to discretize the diffusive term.');
                                break;
                            else
                                df  = 1;
                                Inv = inv(Df);
                                Tf  = df*Inv;
                                cf  = Tf(n,:);
                            end
                        else
                            df      = zeros(1,length(b.xt));
                            df(1,i) = 1;
                            Inv     = inv(Df);
                            cf      = df*Inv;
                        end
                        %  > Update 'x' field.
                        x.cf{a,i}   = cf;
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
                        m.Af{i}(a,b) = x.cf{a,i}(1:length(b));
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
        function [m] = Update_mB(msh,pde,m,s,u,x)
            % >> Set face equation(s) and construct/update B (face dependent coefficient matrix w/ source term contribution).
            for i = 1:size(u.s,2)
                %  > Bf (face contributions).
                if ~isempty(u.s{i})
                    for j = 1:size(u.s{i},1)
                        %  > Face index.
                        a = u.s{i}(j);
                        if ~isempty(s.bv{a,i})
                            %  > Overwrite/update field...
                            m.Bf{i}(j,1) = x.cf{a,i}(end).*s.bv{a,i};
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
                    m.Bt(e,1) = pde.fn.vol(e);
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
        %  > Update 'obj.x.xc' field (nodal solution).
        function [c] = Update_xc(At,Bt)
            c = At\Bt;
        end
        % >> 4.2. ---------------------------------------------------------
        %  > Update 'obj.x.xv' field (nodal values used to reconstruct face).
        function [x] = Update_xv(s,u,x)
            for i = 1:size(u.s,2)
                if ~isempty(u.s{i})
                    for j = 1:size(u.s{i},1)
                        k = u.s{i}(j);
                        %  > Overwrite/update field...
                        if ~isempty(s.c{k,i})
                            %  > Cell contribution(s).
                            l   = s.c{k,i};
                            v.a = x.nv.a.c(l,1);
                            v.x = x.nv.x.c(l,1);
                            %  > Boundary contribution(s).
                            if ~isempty(s.bv{k,i})
                                v.a = [v.a;s.bv{k,i}];
                                v.x = [v.x;s.bv{k,i}];
                            end
                        else
                            v.a = s.bv{k,i};
                            v.x = s.bv{k,i};
                        end
                        x.vf.a{k,i} = v.a;
                        x.vf.x{k,i} = v.x;
                    end
                end
            end
        end
        % >> 4.3. ---------------------------------------------------------
        %  > Update 'pde.x.f' field (reconstructed face values).
        function [x] = Update_xf(u,x)
            for i = 1:size(u.s,2)
                if ~isempty(u.s{i})
                    for j = 1:size(u.s{i},1)
                        k           = u.s{i}(j);
                        x.xf.a(k,i) = x.cf{k,i}*x.vf.a{k,i};
                        x.xf.x(k,i) = x.cf{k,i}*x.vf.x{k,i};
                    end
                end
            end
        end
        
        %% > 5. -----------------------------------------------------------
        % >> 5.1. ---------------------------------------------------------
        %  > Update 'pde.e.p{i}.t.(...)' field (predicted/estimated cell/face truncation error distribution/norms).
        function [e_pt] = Update_et_p(msh,e_pt,sp_l,xp_l,xp_r)
            %  > Error distribution.
            [m,n] = size(e_pt.f);
            for j = 1:n-1
                e_pt.f  (:,j) = sp_l.v(j).*(xp_r.xf.x(:,j)-xp_l.xf.x(:,j));
            end
            e_pt.f      (:,n) = sum(e_pt.f(:,1:n-1),2);
            e_pt.f_abs        = abs(e_pt.f);
            a                 = 1:m-1;
            b                 = a+1;
            e_pt.c      (:,1) = e_pt.f(a,n)-e_pt.f(b,n);
            %  > ... divide cell truncation error by cell volume.
            e_pt.c      (:,1) = e_pt.c;
            e_pt.c_abs  (:,1) = abs(e_pt.c(:,1));
            %  > Error norms.
            e_pt.n.f          = Tools_1D.n(e_pt.f);
            e_pt.n_abs.f      = Tools_1D.n(e_pt.f_abs);
            e_pt.n.c          = Tools_1D.n(e_pt.c,msh.c.Vc);
            e_pt.n_abs.c      = Tools_1D.n(e_pt.c_abs,msh.c.Vc);
        end
        % >> 5.2. ---------------------------------------------------------
        %  > Update 'pde.e.p{i}.c.(...)' field (predicted/estimated cell global discretization error distribution/norms).
        function [e_pc] = Update_ec_p(msh,e_pc,e_pt,m)
            %  > Error distribution.
            e_pc.c    (:,1) = inv(m.At)*e_pt.c;
            %  > ... multiply cell global discretization error by cell volume.
            e_pc.c_abs(:,1) = abs(e_pc.c(:,1));
            %  > Error norms.
            e_pc.n          = Tools_1D.n(e_pc.c,msh.c.Vc);
            e_pc.n_abs      = Tools_1D.n(e_pc.c_abs,msh.c.Vc);
        end
        % >> 5.3. ---------------------------------------------------------
        %  > Update 'pde.e.a{i}.t.(...)' field (analytic cell/face truncation error distribution/norms).
        function [e_at] = Update_et_a(msh,e_at,s,u,x)
            %  > Error distribution.
            [~,n] = size(u.s);
            for i = 1:n
                if ~isempty(u.s{i})
                    for j = 1:size(u.s{i},1)
                        k           = u.s{i}(j);
                        e_at.f(k,i) = s.v(i).*(x.nv.a.f(k,i)-x.xf.a(k,i));
                    end
                    l{i}(:,1) = unique([msh.f.c{[u.s{i}]}]);
                end
            end
            a                   = unique(cat(1,u.s{:}));
            b                   = unique(cat(1,l{:}));
            c                   = b+1;
            e_at.f      (a,n+1) = sum(e_at.f(a,1:n),2);
            e_at.f_abs  (a,:)   = abs(e_at.f(a,:));
            e_at.c      (b,1)   = e_at.f(b,n+1)-e_at.f(c,n+1);
            %  > ... divide cell truncation error by cell volume.
            e_at.c      (b,1)   = e_at.c(b,1);
            e_at.c_abs  (b,1)   = abs(e_at.c(b,1));
            %  > Error norms.
            e_at.n.f            = Tools_1D.n(e_at.f);
            e_at.n_abs.f        = Tools_1D.n(e_at.f_abs);
            e_at.n.c            = Tools_1D.n(e_at.c,msh.c.Vc);
            e_at.n_abs.c        = Tools_1D.n(e_at.c_abs,msh.c.Vc);
        end
        % >> 5.4. ---------------------------------------------------------
        %  > Update 'pde.e.a{i}.c.(...)' field (analytic cell global discretization error distribution/norms).
        function [e_ac] = Update_ec_a(msh,e_ac,x)
            %  > Error distribution.
            e_ac.c    (:,1) = x.nv.a.c-x.nv.x.c;
            e_ac.c_abs(:,1) = abs(e_ac.c);
            %  > Error norms.
            e_ac.n          = Tools_1D.n(e_ac.c,msh.c.Vc);
            e_ac.n_abs      = Tools_1D.n(e_ac.c_abs,msh.c.Vc);
        end
        % >> 5.5. ---------------------------------------------------------
        %  > Update 'pde.e.d{i}.t.(...)' field.
        function [e_dt] = Update_et_d(el_at,er_at,e_dt,e_pt)
            %  > Error distribution.
            n            = 1:size(e_dt.f,2);
            da_f         = el_at.f(:,n)-er_at.f(:,n);
            da_f_abs     = abs(da_f);
            e_dt.f       = e_pt.f(:,n)-da_f;
            e_dt.f_abs   = abs(e_dt.f);
            %  > Error norms.
            e_dt.n.f     = Tools_1D.n(e_dt.f);
            e_dt.n_abs.f = Tools_1D.n(e_dt.f_abs);           
        end
        % >> 5.6 ---------------------------------------------------------
        %  > Update 'pde.e.(...)' field.
        function [e] = Update_e(inp,msh,pde,e,m,s,u,x,add)
            %  > Auxiliary variables.
            for i = 1:size(u.s,2)
                all_s{i}(:,1) = 1:msh.f.Nf;
            end
            ns = inp.pa.ns;
            
            % >> Assign/update fields 'm', 's', 'u' and 'x'.
            i     = 1;
            sp{i} = s;
            mp{i} = m;
            up{i} = u;
            xp{i} = x;
            for i = 1:ns
                up{i+1}        = B_2_1_1D.Set_upd_p(up{i},all_s);
                sp{i+1}        = B_2_1_1D.Update_1 (inp,msh,pde,sp{i},up{i+1},add);
                xp{i+1}        = B_2_1_1D.Update_2 (inp,msh,sp{i+1},up{i+1},xp{i});
                mp{i+1}        = B_2_1_1D.Update_3 (msh,pde,mp{i},sp{i+1},up{i+1},xp{i+1});
                xp{i+1}.nv.x.c = x.nv.x.c;
                xp{i+1}        = B_2_1_1D.Update_4 (sp{i+1},up{i+1},xp{i+1});
            end
            % >> pde.e.p{i}.
            %  > pde.e.p{i}.t(...).
            %  > pde.e.p{i}.c(...).
            for i = 1:ns
                e.p{i}.t = B_2_1_1D.Update_et_p(msh,e.p{i}.t,sp{i},xp{i},xp{i+1});
                e.p{i}.c = B_2_1_1D.Update_ec_p(msh,e.p{i}.c,e.p{i}.t,mp{i});
            end
            % >> pde.e.a{i}.
            if ~inp.pa.comp_av
                nf = 1;
            else
                nf = ns;
            end
            %  > pde.e.a{i}.t(...).
            %  > pde.e.a{i}.c(...).
            for i = 1:nf+1
                if inp.pa.comp_av && i ~= 1
                    xp{i}.nv.x.c = B_2_1_1D.Update_xc(mp{i}.At,mp{i}.Bt);
                end
                e.a{i}.t = B_2_1_1D.Update_et_a(msh,e.a{i}.t,sp{i},up{i},xp{i});
                e.a{i}.c = B_2_1_1D.Update_ec_a(msh,e.a{i}.c,xp{i});
            end  
            % >> pde.e.d{i}.
            %  > pde.e.d{i}.t(...).
            for i = 1:nf
                e.d{i}.t = B_2_1_1D.Update_et_d(e.a{i}.t,e.a{i+1}.t,e.d{i}.t,e.p{i}.t);
            end
        end
        %  > 5.6.1. -------------------------------------------------------
        %  > Auxiliary function (increase method's order).
        function [u] = Set_upd_p(u,u_s)
            A = 2;
            for i = 1:size(u_s,2)
                if ~isempty(u_s{i})
                    for j = 1:size(u_s{i},1)
                        k        = u_s{i}(j);
                        l        = i*size(u_s,2)-1;
                        u.p(k,l) = u.p(k,l)+A;
                    end
                end
                u.s{i} = u_s{i};
            end
        end
        
        %% > 6. -----------------------------------------------------------
        %  > Update/set structure fields.
        function [obj,msh] = Set_struct(obj,msh)
            obj = Tools_1D.Sort_obj(obj);
            msh = Tools_1D.Sort_msh(msh);
        end
    end
end
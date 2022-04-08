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
        %  > Update nodal/face values.
        function [x] = Update_4(s,u,x,fld_u)
            x        = B_2_1_1D.Update_xv(s,u,x,fld_u);
            x        = B_2_1_1D.Update_xf(u,x,fld_u);
        end
        % >> 1.5. ---------------------------------------------------------
        %  > Update '1', '2' and '3'.
        function [s,x,m] = Update_123(inp,msh,f,m,s,u,x,add_b)
            s            = B_2_1_1D.Update_1(inp,msh,f,s,u,add_b);
            x            = B_2_1_1D.Update_2(msh,s,u,x);
            m            = B_2_1_1D.Update_3(msh,f,m,s,u,x);
        end        
        % >> 1.6. ---------------------------------------------------------
        %  > Update all...
        function [e,f,s,u,x] = Update_all(inp,msh,e,f,m,s,u,x,add_b,fld_u)
            % >> Update fields 'm', 's' and 'x'.
            [s,x,m] = B_2_1_1D.Update_123(inp,msh,f,m,s,u,x,add_b);
            % >> Update field 'x'.
            x.nv.x.c = B_2_1_1D.Update_xc(m.At,m.Bt);
            x        = B_2_1_1D.Update_4 (s,u,x,fld_u);
            % >> Update field 'e'.
            e        = B_2_1_1D.Update_e (inp,msh,e,f,m,s,u,x,add_b);
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
        function [x] = Update_xv(s,u,x,fld_u)
            for i = 1:length(fld_u)  
                j = fld_u(i);
                for k = 1:size(u.s,2)
                    if ~isempty(u.s{k})
                        for l = 1:size(u.s{k},1)
                            m = u.s{k}(l);
                            %  > Overwrite/update field...
                            if ~isempty(s.c{m,k})
                                %  > Cell contribution(s).
                                n     = s.c{m,k};
                                v.(j) = x.nv.(j).c(n,1);
                                %  > Boundary contribution(s).
                                if ~isempty(s.bv{m,k})
                                    v.(j) = [v.(j);s.bv{m,k}];
                                end
                            else
                                v.(j) = s.bv{m,k};
                            end
                            x.vf.(j){m,k} = v.(j);
                        end
                    end
                end
            end
        end
        % >> 4.3. ---------------------------------------------------------
        %  > Update 'x.cf', 'x.ff' and 'x.xf' (fitted face coefficients/function handles/values).
        function [x] = Update_xf(u,x,fld_u)
            for i = 1:length(fld_u)
                j = fld_u(i);
                for k = 1:size(u.s,2)
                    if ~isempty(u.s{k})
                        for l = 1:size(u.s{k},1)
                            m = u.s{k}(l);
                            %  > 'x.xf' and 'x.cf'.
                            x.xf.(j)(m,k) = x.Tf{m,k}*x.vf.(j){m,k};
                            x.cf.(j){m,k} = x.if{m,k}*x.vf.(j){m,k};
                            %  > 'x.ff'.
                            %x.ff.(j){m,k} = B_2_1_1D.fit(x.cf.(j){m,k},k);  
                        end
                    end
                end
            end
        end

        %% > 5. -----------------------------------------------------------
        % >> 5.1. ---------------------------------------------------------
        %  > Update 'e.p(...)' field (predicted/estimated cell/face truncation error distribution/norms).
        function [ep] = Update_ep(ep,mp_o,sp_o,xp_o,xp_n,Vc)
            %  > \tau_f(\phi) & \tau(\nabla\phi).
            for i = 1:size(ep.t.f,2)-1
                ep.t.f(:,i) = sp_o.v(i).*(xp_n.xf.x(:,i)-xp_o.xf.x(:,i));
            end
            %  > Update remaining error fields...
            ep = Tools_1D.Set_e(ep,mp_o,Vc); 
        end
        % >> 5.2. ---------------------------------------------------------
        %  > Update 'e.a(...)' field (analytic cell/face truncation error distribution/norms).
        function [ea] = Update_ea(ea,mp_o,sp_o,xp_o,Vc)
            %  > \tau_f(\phi) & \tau(\nabla\phi).
            for i = 1:size(ea.t.f,2)-1
                ea.t.f(:,i) = sp_o.v(i).*(xp_o.nv.a.f(:,i)-xp_o.xf.a(:,i));
            end
            %  > Update remaining error fields...
            ea = Tools_1D.Set_e(ea,mp_o,Vc); 
        end
        % >> 5.3. ---------------------------------------------------------
        %  > Update 'e.d(...)' field.
        function [ed] = Update_ed(ea,ed,ep,Vc)
            %  > Field.
            f         = ["c","t"];
            g{1}(1)   = "c";
            g{2}(1,:) = ["c","f"];
            %  > Error distribution.
            for i = 1:size(g,2)
                a = f(i);
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
            for i = 1:size(u.s,2)
                all_s{i}(:,1) = 1:msh.f.Nf;
            end
            fld = ["a","x"];

            for i = 1:inp.pa.ns
                %  > Initialize...
                if i == 1
                    mp_o = m;
                    sp_o = s;
                    up_o = u;
                    xp_o = x;
                    ea_o = e.a{i};
                    ed_o = e.d{i};
                    ep_o = e.p{i};
                    ea_n = ea_o;
                    ed_n = ed_o;
                    ep_n = ep_o;
                end
                % >> Assign/update fields 'e', 'm', 's', 'u' and 'x'.
                %  > 'm', 's', 'u' and 'x'.
                
                
                
                
                
                
                up_n = B_2_1_1D.Set_upd_p(up_o,all_s);
                sp_n = B_2_1_1D.Update_1 (inp,msh,f,sp_o,up_n,add_b);
                xp_n = B_2_1_1D.Update_2 (msh,sp_n,up_n,xp_o);
                mp_n = B_2_1_1D.Update_3 (msh,f,mp_o,sp_n,up_n,xp_n);
                xp_n = B_2_1_1D.Update_4 (sp_n,up_n,xp_n,fld); 
                
                %  > e.
                Vc     = msh.c.Vc;
                ea_n   = B_2_1_1D.Update_ea(ea_n,mp_o,sp_o,xp_o,Vc);
                ep_n   = B_2_1_1D.Update_ep(ep_n,mp_o,sp_o,xp_o,xp_n,Vc);
                ed_n   = B_2_1_1D.Update_ed(ea_n,ed_n,ep_n,Vc);
                % >> Assign...
                %  > Update next cycle.
                if i ~= inp.pa.ns
                    mp_o = mp_n;
                    sp_o = sp_n;
                    up_o = up_n;
                    xp_o = xp_n;
                    ea_o = ea_n;
                    ed_o = ed_n;
                    ep_o = ep_n;
                end
                %  > e.
                e.a{i} = ea_n;
                e.d{i} = ed_n;
                e.p{i} = ep_n;
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
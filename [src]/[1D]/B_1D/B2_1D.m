classdef B2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Update matrices.
        function [m] = Update_3(inp,msh,f,m,s,u,x)
            m        = B2_1D.Update_mA(inp,msh,f,m,s,u,x);
            m        = B2_1D.Update_mB(inp,msh,f,m,s,u,x);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update nodal/face values, etc.
        function [x] = Update_4(f,s,u,x)
            x        = B2_1D.Update_xv(f,s,u,x);
            x        = B2_1D.Update_cf(u,x);
            x        = B2_1D.Update_xf(u,x);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Update matrix A.
        function [m] = Update_mA(inp,msh,f,m,s,u,x)
            %  > Auxiliary variables.
            Nc = msh.c.Nc;
            
            % >> Set face equation(s) and update A (cell dependent coefficient matrix).
            for i = 1:size(u.s,2)
                %  > Af (face contributions).
                if ~isempty(u.s{i})
                    for j = u.s{i}'
                        k                   = ismembc(s.t{j,i},f.bd.x);
                        m.Af{i}(j,:)        = zeros(1,Nc);
                        m.Af{i}(j,s.c{j,i}) = inp.c(i).*x.Tf{j,i}(~k);
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
                m.nnz.Ac(i) = B2_1D.Set_nnz(m.Ac{i});
                m.nnz.Af(i) = B2_1D.Set_nnz(m.Af{i});
            end
            %  > Cell(s) to be updated.
            b = unique(cat(1,a{:}));
            
            %  > At (cumulative/total matrix).
            if ~isempty(b)
                for i = b'
                    %  > Add convective/diffusive contributions...
                    m.At(i,:) = zeros(1,Nc);
                    for j = 1:size(u.s,2)
                        m.At(i,:) = m.At(i,:)+m.Ac{j}(i,:);
                    end
                end
            end
            m.nnz.At = B2_1D.Set_nnz(m.At);
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Update matrix B.
        function [m] = Update_mB(inp,msh,f,m,s,u,x)
            %  > Auxiliary variables.
            Nf   = msh.f.Nf;
            k_bf = zeros(Nf,size(u.s,2));
            
            % >> Set face equation(s) and update B (face dependent coefficient matrix w/ source term contribution).
            for i = 1:size(u.s,2)
                %  > Bf (face contributions).
                if ~isempty(u.s{i})
                    for j = u.s{i}'
                        k = ismembc(s.t{j,i},f.bd.x);
                        if any(k)
                            k_bf   (j,i) = 1;
                            m.Bf{i}(j,1) = inp.c(i).*x.Tf{j,i}(k).*f.bd.v(f.bd.x == s.t{j,i}(k));
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
                    m.Bt(i,1) = f.st(i);
                    for j = 1:size(u.s,2)
                        m.Bt(i,1) = m.Bt(i,1)+m.Bc{j}(i,1);
                    end
                end
            end
        end
        % >> 2.3. ---------------------------------------------------------
        %  > Update 'nnz' of each matrix.
        function [nnz_m] = Set_nnz(m)
            nnz_m = nnz(m);
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > Update 'x.nv.x.c' field (nodal solution).
        function [xc] = Update_xc(At,Bt)
            xc = At\(Bt);
        end
        % >> 3.2. ---------------------------------------------------------
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
        % >> 3.3. ---------------------------------------------------------
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
        % >> 3.4. ---------------------------------------------------------
        %  > Update 'x.f' field (face values).
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
        
        %% > 4. -----------------------------------------------------------
        % >> 4.1. ---------------------------------------------------------
        %  > Update 'e.p(...)' field (predicted/estimated cell/face truncation error distribution/norms).
        function [ed,ep,x] = Update_ed_ep(inp,eda,ed,ep,f,m,s,u,x,Vc,fs)
            %  > \tau_f(\phi), \tau(\nabla\phi), \tau_f and \tau_c.
            ep = B2_1D.Set_1_e(ep,inp.c,x{1}.xf.x,x{2}.xf.x);
            %  > Update remaining error fields...
            if ~fs(2)
                if ~fs(1)
                    %  > w/ LO (j=2).
                    j = 2;
                else
                    %  > w/ HO (j=1).
                    j = 1;
                end
                ep = B2_1D.Set_2_e(ep,m{j},Vc);
            else
                %  > Add as source(?).
                j  = 1;
                x  = B2_1D.Add_ep_tc(ep,f,m{j},s,u,x);
                ep = B2_1D.Set_1_e  (ep,s{1}.v,x{1}.xf.x,x{2}.xf.x);
                ep = B2_1D.Set_2_e  (ep,m{j},Vc);
            end
            ed = B2_1D.Update_ed(eda,ed,ep,Vc);
        end
        %  > 4.1.1. -------------------------------------------------------
        %  > Auxiliary function (add \tau_c(p) as source to the RHS).
        function [x] = Add_ep_tc(ep,f,m,s,u,x)
            %  > Update cell values...
            xc = B2_1D.Update_xc(m.At,m.Bt+ep.t.c);
            %  > Update "x.(...).x" field only.
            for i = 1:size(x,2)
                u{i}.f      = "x";
                x{i}.nv.x.c = xc;
                x{i}        = B2_1D.Update_4(f,s{i},u{i},x{i});
            end
        end
        % >> 4.2. ---------------------------------------------------------
        %  > Update 'e.a(...)' field (analytic cell/face truncation error distribution/norms).
        function [ea] = Update_ea(inp,msh,ea,f,m,s,u,x,Vc)
            %  > \tau_f(\phi), \tau(\nabla\phi), \tau_f and \tau_c.
            ea = B2_1D.Set_1_e(ea,inp.c,x.xf.a,x.nv.a.f);
            %  > Update remaining error fields...
            ea = B2_1D.Set_2_e(ea,m,Vc);
            %  > Update truncated terms' magnitude (w/ anlytic field).
            if inp.t_terms.allow && ~inp.p_adapt.allow
                nt   = inp.t_terms.n;
                ea.m = B1_1D.Update_t_terms(inp,msh,ea.m,f,nt,s,u,x);
            end 
        end
        % >> 4.3. ---------------------------------------------------------
        %  > Update 'e.da(...)' field (analytic cell/face (difference) truncation error distribution/norms).
        function [eda] = Update_eda(inp,eda,m,x,Vc)
            %  > \tau_f(\phi), \tau(\nabla\phi), \tau_f and \tau_c.
            eda = B2_1D.Set_1_e(eda,inp.c,x{1}.xf.a,x{2}.xf.a);
            %  > Update remaining error fields...
            eda = B2_1D.Set_2_e(eda,m,Vc);
        end
        % >> 4.4. ---------------------------------------------------------
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
            ed.c.n       = B2_1D.Set_n(ed.c.c,Vc);
            ed.c.n_abs   = B2_1D.Set_n(ed.c.c_abs,Vc);
            ed.t.n.c     = B2_1D.Set_n(ed.t.c,Vc);
            ed.t.n.f     = B2_1D.Set_n(ed.t.f);
            ed.t.n_abs.c = B2_1D.Set_n(ed.t.c_abs,Vc);
            ed.t.n_abs.f = B2_1D.Set_n(ed.t.f_abs);
        end
        % >> 4.5. ---------------------------------------------------------
        %  > Auxiliary functions.
        %  > 4.5.1. -------------------------------------------------------
        %  > Auxiliary function #1.
        %  > Compute error norms (cell/face L1,L2 and L_infinity norms).
        function [L] = Set_n(E,V)
            if nargin == 1
                L(1,:) = mean(E);
                L(2,:) = mean(sqrt(E.^2));
                L(3,:) = max(E);
            else
                L(1,:) = sum (E.*V)./sum(V);
                L(2,:) = sum (sqrt((E.*V).^2))./sum(sqrt(V.^2));
                L(3,:) = max(E);
            end
        end
        %  > 4.5.2. -------------------------------------------------------
        %  > Auxiliary function #2.
        %  > Compute remaining error fields (based on the convective/diffusive facial components).
        function [e] = Set_1_e(e,v,x,y)
            %  > \tau_f(\phi) and \tau(\nabla\phi).
            [a,b] = size(e.t.f);
            for i = 1:b-1
                e.t.f(:,i) = v(i).*(y(:,i)-x(:,i));
            end
            %  > \tau_f.
            e.t.f    (:,b) = sum(e.t.f(:,1:b-1),2);
            %  > \tau_c.
            for i = 1:a-1
                e.t.c(i,1) = e.t.f(i,b)-e.t.f(i+1,b);
            end
        end
        %  > 4.5.3. -------------------------------------------------------
        %  > Auxiliary function #3.
        function [e] = Set_2_e(e,m,Vc)
            % >> Error distribution.
            %  > e_c.
            e.c.c    (:,1) = m.At\e.t.c;
            %  > abs().
            e.c.c_abs      = abs(e.c.c);
            e.t.c_abs      = abs(e.t.c);
            e.t.f_abs      = abs(e.t.f);
            % >> Error norms.
            e.c.n          = B2_1D.Set_n(e.c.c,Vc);
            e.c.n_abs      = B2_1D.Set_n(e.c.c_abs,Vc);
            e.t.n.f        = B2_1D.Set_n(e.t.f);
            e.t.n_abs.f    = B2_1D.Set_n(e.t.f_abs);
            e.t.n.c        = B2_1D.Set_n(e.t.c,Vc);
            e.t.n_abs.c    = B2_1D.Set_n(e.t.c_abs,Vc);
        end 
        % >> 4.6 ----------------------------------------------------------
        %  > Update 'e.(...)' field (error).
        function [e,m,s,x] = Update_e(inp,msh,e,f,m,s,u,x,f_xc)
            %  > Auxiliary variables.
            nc    = size(x,2);
            Nf    = msh.f.Nf;
            Vc    = msh.c.Vc;
            ch    = inp.b.change   (2); %  > Boundary treatment.
            fs(1) = inp.p_adapt.opt(1); %  > Use higher-order solution(?).
            fs(2) = inp.p_adapt.opt(2); %  > Add lower-order (predicted) cell truncation error as source term to the RHS(?).
            
            % >> #1: Update stencil/nodal solution...
            %  > Update fields 'm', 's' and 'x'.
            ic                  = 2;
            [m{ic},s{ic},x{ic}] = B3_1D.Update_all(inp,msh,f,m{ic},s{ic},u{ic},x{ic},ch,[1,fs(1),fs(1)]);
            %  > Update nodal solution...
            if f_xc
                if ~fs(1)
                    ic = 1;
                    xc = B2_1D.Update_xc(m{ic}.At,m{ic}.Bt);
                    for i = 1:size(x,2)
                        x{i}.nv.x.c = xc;
                    end
                else
                    ic           = 1;
                    x{ic}.nv.x.c = x{ic+1}.nv.x.c;
                end
            end
            %  > Update field 'x'.
            for i = 1:size(x,2)
                x{i} = B2_1D.Update_4(f,s{i},u{i},x{i});
            end
            % >> #2: Update fields 'e.a' and 'e.da': if the (predicted) cell truncation error is added as a source term, no need to update matrices, since e=A(LO)\(\tau_c).
            for i = 1:nc
                e.a{i} = B2_1D.Update_ea(inp,msh,e.a{i},f,m{i},s{i},u{i},x{i},Vc);
                if i ~= 1
                    j         = i-1:i;
                    e.da{i-1} = B2_1D.Update_eda(inp,e.da{i-1},m{i},x(j),Vc);
                end
            end
            % >> #3: Update fields 'e.d', 'e.p' and 'x'.
            i                    = 1;
            j                    = i:i+1;
            [e.d{i},e.p{i},x(j)] = ...
                B2_1D.Update_ed_ep(inp,e.da{i},e.d{i},e.p{i},f,m(j),s(j),u(j),x(j),Vc,fs);
        end
    end
end
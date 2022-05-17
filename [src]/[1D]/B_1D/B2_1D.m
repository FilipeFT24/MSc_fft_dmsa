classdef B2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Update matrices.
        function [m] = Update_3(inp,msh,f,m,s,u,x)
            m        = B2_1D.Update_m(inp,msh,f,m,s,u,x);
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
        %  > Update matrices A and B.
        function [m] = Update_m(inp,msh,f,m,s,u,x)
            %  > Initialize rows (r) to be updated...
            uf        = [msh.f.ic{RunLength(sort(cat(1,u.s{:})))}]';
            r         = 1:100;%RunLength(uf)';
            m.At(r,:) = 0;
            m.Bt(r)   = f.st(r);
            for i = 1:numel(u.s)
                m.Ac{i}(r,:) = 0;
                m.Bc{i}(r)   = 0;
            end
            
            %  > For each term (convective/diffusive)...
            for i = 1:numel(u.s)
                for j = r
                    for k = 1:numel(msh.c.f.if(j,:)), ff = msh.c.f.if(j,k);
                        %  > Cell/face indices used to fit face "ff"...
                        l = s.logical{ff,i};
                        %  > Auxiliary variables.
                        a = s.i{ff,i}( l);
                        b = s.i{ff,i}(~l);
                        switch k
                            case 1
                                Nf = -1;
                            case 2
                                Nf =  1;
                            otherwise
                                return;
                        end
                        %  > Ac (cell contributions).
                        m.Ac{i}(j,a) = m.Ac{i}(j,a)+Nf*x.Tf_V{ff,i}(l);
                        %  > Bc (cell contributions).
                        if any(~l)
                            m.Bc{i}(j,1) = m.Bc{i}(j)-Nf*x.Tf_V{ff,i}(~l)*f.bd.v(ismembc(f.bd.i,sort(b)));
                        end
                    end
                end
                %  > At and Bt (cumulative/total matrices).
                m.At = m.At+m.Ac{i};
                m.Bt = m.Bt+m.Bc{i};
            end
            %  > nnz.
            for i = 1:numel(u.s)
                m.nnz.Ac(i) = nnz(m.Ac{i});
            end
            m.nnz.At = nnz(m.At);
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
                        for k = 1:101%u.s{j}'
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
                        for k = 1:101%u.s{j}'
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
                        for k = 1:101%u.s{j}'
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
            ep = B2_1D.Set_1_e(ep,inp.c,x{1}.xf.x,x{2}.xf.x,x{2}.xf.x);
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
                ep = B2_1D.Set_1_e  (ep,s{1}.v,x{1}.xf.x,x{2}.xf.x,x{2}.xf.x);
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
            ea = B2_1D.Set_1_e(ea,m,inp.c,x.xf.a,x.xf.x); %x.nv.a.f
%             
%             ea.c.c    (:,1) = m.At\ea.t.c;
            ea.c.c = x.nv.x.c-f.av.c;
%             disp(mean(abs(ea.c.c-sa)));
%             
%             
%             figure(1);
%             hold on;
%             plot(abs(ea.c.c),'-r');
%             plot(abs(sa),'ob');
%             set(gca,'YScale','log');
%             
            
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
            eda = B2_1D.Set_1_e(eda,inp.c,x{1}.xf.a,x{2}.xf.a,x{2}.xf.a);
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
        function [e] = Set_1_e(e,m,v,x,y)
            %  > \tau_f(\phi) and \tau(\nabla\phi).
            [a,b] = size(e.t.f);
            for i = 1:b-1
                e.t.f(:,i) = v(i).*(x(:,i)-y(:,i));
            end
            %  > \tau_f.
            e.t.f(:,b) = sum(e.t.f(:,1:b-1),2);
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
            %ic                  = 2;
            %[m{ic},s{ic},x{ic}] = B3_1D.Update_all(inp,msh,f,m{ic},s{ic},u{ic},x{ic},ch,[1,fs(1),fs(1)]);
            %  > Update nodal solution...
            if f_xc
                if ~fs(1)
                    ic = 1;
                    xc = B2_1D.Update_xc(m{ic}.At,m{ic}.Bt);
                    for i = 1:1%size(x,2)
                        x{i}.nv.x.c = xc;
                    end
                else
                    ic           = 1;
                    x{ic}.nv.x.c = x{ic+1}.nv.x.c;
                end
            end
            %  > Update field 'x'.
            for i = 1:1%size(x,2)
                x{i} = B2_1D.Update_4(f,s{i},u{i},x{i});
            end
            % >> #2: Update fields 'e.a' and 'e.da': if the (predicted) cell truncation error is added as a source term, no need to update matrices, since e=A(LO)\(\tau_c).
            for i = 1:1%nc
                e.a{i} = B2_1D.Update_ea(inp,msh,e.a{i},f,m{i},s{i},u{i},x{i},Vc);
%                 if i ~= 1
%                     j         = i-1:i;
%                     e.da{i-1} = B2_1D.Update_eda(inp,e.da{i-1},m{i},x(j),Vc);
%                 end
            end
            e.a{1}.c.c_abs = mean(abs(f.av.c-x{1}.nv.x.c));
            %e.a{1}.c.n_abs = B2_1D.Set_n(e.a{1}.c.c_abs,Vc);
            
            
%             % >> #3: Update fields 'e.d', 'e.p' and 'x'.
%             i                    = 1;
%             j                    = i:i+1;
%             [e.d{i},e.p{i},x(j)] = ...
%                 B2_1D.Update_ed_ep(inp,e.da{i},e.d{i},e.p{i},f,m(j),s(j),u(j),x(j),Vc,fs);
         end
    end
end
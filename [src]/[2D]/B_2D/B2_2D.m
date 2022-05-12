classdef B2_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "m" (matrices).
        function [m] = Initialize_3(nc,ns,Nc)
            for i = 1:ns
                for j = 1:nc(2)
                    m{i}.Ac{j} = zeros(Nc);
                    m{i}.Bc{j} = zeros(Nc);
                end
                m{i}.At = zeros(Nc);
                m{i}.Bt = zeros(Nc,1);
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update matrices.
        function [m] = Update_3(inp,msh,f,m,s,u,x)
            m        = B2_2D.Update_m(inp,msh,f,m,s,u,x);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Update matrices A and B.
        function [m] = Update_m(inp,msh,f,m,s,u,x)
            %  > Initialize rows (r) to be updated...
            for i = 1:numel(u.s)
                uf{i} = [msh.f.ic{RunLength(sort(cat(1,u.s{i}{:})))}];
            end
            r = RunLength(sort(cat(2,uf{:})));
            for i = 1:numel(u.s)
                m.Ac{i}(r,:) = 0;
                m.Bc{i}(r)   = 0;
            end
            m.At(r,:) = 0;
            m.Bt(r)   = f.st(r);

            %  > For each term (convective/diffusive)...
            for i = 1:numel(u.s)
                for j = r
                    for k = 1:numel(msh.c.f.if(j,:)), ff = msh.c.f.if(j,k);
                        %  > Cell/face indices used to fit face "ff" (for each direction)...
                        l = s.logical{ff,i};
                        for n = 1:numel(l)
                            %  > Auxiliary variables.
                            a            = s.i{ff,i} {n}( l{n});
                            b            = s.i{ff,i} {n}(~l{n});
                            Sf_n         = msh.c.f.Sf{j}( k,n);
                            %  > Ac (cell contributions).
                            m.Ac{i}(j,a) = m.Ac{i}(j,a)+Sf_n*x.Tf_V{ff,i}{n}(:,l{n});
                            %  > Bc (cell contributions).
                            if any(~l{n})
                                m.Bc{i}(j) = m.Bc{i}(j)-Sf_n*x.Tf_V{ff,i}{n}(:,~l{n})*f.bd.v(ismembc(f.bd.i,sort(b)));
                            end
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
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Initialize field "e" (error).
        function [e] = Initialize_e(nc,ns,Nc,Nf)
            %  > Fields: 'e.a', 'e.d' and 'e.p'.
            for i = ["a","da","p","d"]
                switch i
                    case "a"
                        l = ns;
                    otherwise
                        l = ns-1;
                end
                for j = 1:l
                    %  > Error distribution.
                    e.(i){j}.c.c_abs   = zeros(Nc,1);
                    e.(i){j}.t.c       = zeros(Nc,1);
                    e.(i){j}.t.c_abs   = zeros(Nc,1);
                    e.(i){j}.t.f_abs   = zeros(Nf,nc(1)+1);
                    %  > Error norms.
                    e.(i){j}.c.n_abs   = zeros(3 ,1);
                    e.(i){j}.t.n_abs.c = zeros(3 ,1);
                    e.(i){j}.t.n_abs.f = zeros(3 ,nc(1)+1);
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update 'e.(...)' field (error).
        function [e] = Update_e(inp,msh,e,f,m,s,u,x)
            %  > Auxiliary variables.
            nc = 1;
            Vc = msh.c.Volume;
            
            % >> #2: Update fields 'e.a' and 'e.da': if the (predicted) cell truncation error is added as a source term, no need to update matrices, since e=A(LO)\(\tau_c).
            for i = 1:nc
                e.a{i} = B2_2D.Update_ea(inp,msh,e.a{i},f,m{i},s{i},u{i},x{i},Vc);
            end
        end
        %  > 1.2.1. ---------------------------------------------------------
        %  > Update 'e.a(...)' field (analytic cell/face truncation error distribution/norms).
        function [ea] = Update_ea(inp,msh,ea,f,m,s,u,x,Vc)
            %  > \tau_f(\phi), \tau(\nabla\phi), \tau_f and \tau_c.
            ea = B2_2D.Set_1_e(msh,ea,f,x);
            %  > Update remaining error fields...
            ea = B2_2D.Set_2_e(msh,ea,f,x,m,Vc);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Auxiliary functions.
        %  > 1.3.1. -------------------------------------------------------
        %  > Auxiliary function #1.
        %  > Compute remaining error fields (based on the convective/diffusive facial components).
        function [e] = Set_1_e(msh,e,f,x)
            %  > Auxiliary variables.
            n    = size(x.xf.a,2);
            d_tf = cell(1,n+1);
            
            %  > "Average face value".
            %  > \tau_f(\phi)_(x,y) and \tau_f(\nabla\phi)_(x,y).
            for i = 1:n+1
                if i ~= n+1
                    d_tf{i} = x.xf.a{i}-x.xf.x{i};
                else
                    d_tf{i} = sum(cat(3,d_tf{1:n}),3);
                end
            end
            %  > \tau_c and \tau_c_abs.
            for i = 1:msh.c.Nc
                for j = 1:numel(msh.c.f.if(i,:))
                    e.t.c(i) = e.t.c(i)+msh.c.f.Sf{i}(j,:)*d_tf{n+1}(msh.c.f.if(i,j),:)';
                end
            end
            e.t.c_abs = abs(e.t.c);
            %  > \tau_f_abs.
            for i = 1:msh.f.Nf
                k       = 1;
                c       = msh.f.ic  {i}(k);
                Sf(i,:) = msh.c.f.Sf{c}(msh.c.f.if(c,:) == i,:);
                for j = 1:n+1
                    e.t.f_abs(i,j) = abs(d_tf{j}(i,:)*Sf(i,:)');
                end
            end
        end
        %  > 1.3.2. -------------------------------------------------------
        %  > Auxiliary function #2.
        %  > Compute "average face value".
        function [y] = nv_a_f_V(f,gf)
            for i = 1:size(gf,2)
                for j = 1:size(gf,1)
                    for k = 1:numel(gf{j,i})
                        if i == 1
                            y{i}(j,k) = gf{j,i}{k}.V'*f.fh.f.f   (gf{j,i}{k}.xy);
                        else
                            y{i}(j,k) = gf{j,i}{k}.V'*f.fh.f.d{k}(gf{j,i}{k}.xy);
                        end
                    end
                end
            end
        end
        %  > 1.3.3. -------------------------------------------------------
        %  > Auxiliary function #4.
        function [e] = Set_2_e(msh,e,f,x,m,Vc)
            % >> Error distribution.
            %  > e_c_abs.
            e.c.c_abs   = abs(m.At\e.t.c);
            % >> Error norms.
            e.c.n_abs   = Tools_2.Set_n(e.c.c_abs,Vc);
            e.t.n_abs.c = Tools_2.Set_n(e.t.c_abs,Vc);
            e.t.n_abs.f = Tools_2.Set_n(e.t.f_abs);
        end
    end
end
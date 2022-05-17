classdef B2_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "m" (matrices).
        function [m] = Initialize_3(nc,ns,Nc)
            for i = 1:ns
                %  > Matrices.
                for j = 1:nc(2)
                    m{i}.Ac{j} = zeros(Nc);
                    m{i}.Bc{j} = zeros(Nc,1);
                end
                m{i}.At = zeros(Nc);
                m{i}.Bt = zeros(Nc,1);
                %  > nnz.
                m{i}.nnz.Ac = zeros(1,nc(2));
                m{i}.nnz.At = 0;
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
            r         = RunLength(sort(cat(2,uf{:})));
            m.At(r,:) = 0;
            m.Bt(r,1) = f.st(r);
            
            %  > For each term (convective/diffusive)...
            for i = 1:numel(u.s)
                %  > Initialize...
                m.Ac{i}(r,:) = 0;
                m.Bc{i}(r,1) = 0;
                %  > Loop through rows...
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
                                m.Bc{i}(j,1) = m.Bc{i}(j,1)-Sf_n*x.Tf_V{ff,i}{n}(:,~l{n})*f.bd.v(ismembc(f.bd.i,sort(b)));
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
            for i = ["a","p"]
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
                    for k = ["x","y"]
                        e.(i){j}.t.f_abs.  (k) = zeros(Nf,nc(1)+1);
                        e.(i){j}.t.n_abs.f.(k) = zeros(3 ,nc(1)+1);
                    end
                    %  > Error norms.
                    e.(i){j}.c.n_abs   = zeros(3 ,1);
                    e.(i){j}.t.n_abs.c = zeros(3 ,1);
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update 'e.(...)' field (error).
        function [e] = Update_e(inp,msh,e,f,m,s,u,x)
            % >> #2: Update field 'e.a'.
            i      = 1;
            e.a{i} = B2_2D.Update_ea(msh,e.a{i},f,m{i},x{i});
        end
        %  > 1.2.1. ---------------------------------------------------------
        %  > Update 'e.a(...)' field (analytic cell/face truncation error distribution/norms).
        function [ea] = Update_ea(msh,ea,f,m,x)
            %  > Auxiliary variables.
            n      = numel(m.Ac)+1;
            e_t_fx = cell (1,n);
            
            % >> f.
            %  > \tau_f: {1}-\tau_f(\phi).
            %            {2}-\tau_f(\nabla\phi).
            %            {3}-\tau_f.
            for i = 1:n
                if i ~= n
                    e_t_fx{i} = x.xf.a{i}-x.xf.x{i};
                else
                    e_t_fx{i} = sum(cat(3,e_t_fx{1:n-1}),3);
                end
            end
            %  > \tau_f_abs: x(v,g).
            %                y(v,g).
            fn = fieldnames(ea.t.f_abs);
            l  = 1;
            for i = 1:msh.f.Nf
                c       = msh.f.ic  {i}(l);
                Sf(i,:) = msh.c.f.Sf{c}(msh.c.f.if(c,:) == i,:);
                for j = 1:n-1
                    for k = 1:n
                        ea.t.f_abs.(fn{j})(i,k) = abs(e_t_fx{k}(i,j).*Sf(i,j));
                    end
                end
            end
            %  > Error norms.
            for j = 1:n-1
                ea.t.n_abs.f.(fn{j}) = Tools_2.Set_n(ea.t.f_abs.(fn{j}));
            end
            
            % >> c.
            %  > \tau_c.
            for i = 1:msh.c.Nc
                for j = 1:numel(msh.c.f.if(i,:))
                    ea.t.c(i) = ea.t.c(i)+msh.c.f.Sf{i}(j,:)*e_t_fx{n}(msh.c.f.if(i,j),:)';
                end
            end
            %  > \tau_c_abs.
            ea.t.c_abs   = abs(ea.t.c);
            %  > e_c_abs.
            %  > Equivalent to: ea.c.c_abs = f.av.c-x.nv.x.c.
            ea.c.c_abs   = abs(m.At\ea.t.c);
            %  > Error norms.
            ea.c.n_abs   = Tools_2.Set_n(ea.c.c_abs,msh.c.Volume);
            ea.t.n_abs.c = Tools_2.Set_n(ea.t.c_abs,msh.c.Volume);
        end
    end
end
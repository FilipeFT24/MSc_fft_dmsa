classdef B2_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "m" (matrices).
        function [m] = Initialize_3(nc,ns,Nc)
            for i = 1:ns
                for j = 1:nc(2)
                    m{i}.Ac{j} = zeros(Nc);
                    m{i}.Bc{j} = zeros(Nc,1);
                end
                m{i}.At     = zeros(Nc);
                m{i}.Bt     = zeros(Nc,1);
                m{i}.nnz.Ac = zeros(1,nc(2));
                m{i}.nnz.At = 0;
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update matrices.
        function [m] = Update_3(msh,f,m,s,u,x)
            m        = B2_2D.Update_m(msh,f,m,s,u,x);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Update matrices A and B.
        function [m] = Update_m(msh,f,m,s,u,x)
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
                            a    = s.i{ff,i} {n}( l{n});
                            b    = s.i{ff,i} {n}(~l{n});
                            Sf_n = msh.c.f.Sf{j}( k,n);
                            %  > Ac (cell contributions).
                            m.Ac{i}(j,a) = m.Ac{i}(j,a)+Sf_n*x.Tf_V{ff,i}{n}(:,l{n});
                            %  > Bc (cell contributions).
                            if any(~l{n})
                                i_bd         = arrayfun(@(x) find(f.bd.i == x),b);
                                m.Bc{i}(j,1) = m.Bc{i}(j,1)-Sf_n*x.Tf_V{ff,i}{n}(:,~l{n})*f.bd.v(i_bd);
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
        % >> 1.4. ---------------------------------------------------------
        %  > Update nodal solution/face values (multiplied by V).
        function [x] = Update_4(f,m,s,u,x,f_xc)
            %  > Update nodal solution(?).
            if f_xc
                x.nv.x.c = m.At\m.Bt;
            end
            
            %  > Check if empty...
            fn   = convertCharsToStrings(fieldnames(x.vf));
            fn_n = numel(fn);
            if all(cellfun(@isempty,x.vf.a),'all')
                flag = true;
            end
            %  > Update only necessary/all entries for field "a"/"x".
            for i = 1:numel(u.s)
                for j = 1:numel(u.s{i})
                    for k = 1:fn_n
                        switch fn(k)
                            case "a"
                                if ~flag
                                    r = u.s{i}{j}';
                                else
                                    r = 1:size(x.xf.a{i},1);
                                end
                            case "x"
                                r = 1:size(x.xf.a{i},1);
                        end
                        for l = r
                            %  > Cell/face indices used to fit "k".
                            m = s.logical{l,i}{j};
                            a = s.i      {l,i}{j}( m);
                            b = s.i      {l,i}{j}(~m);
                            %  > Cell value(s).
                            x.vf.(fn(k)){l,i}{j}(m,1) = x.nv.(fn(k)).c(s.i{l,i}{j}(m));
                            %  > Face value(s).
                            if any(~m)
                                x.vf.(fn(k)){l,i}{j}(~m,1) = f.bd.v(ismembc(f.bd.i,sort(b)));
                            end
                            %  > "x.cf"
                            x.cf.(fn(k)){l,i}{j}   = x.Pf  {l,i}{j}*x.vf.(fn(k)){l,i}{j};
                            %  > "x.xf".
                            x.xf.(fn(k))  {i}(l,j) = x.Tf_V{l,i}{j}*x.vf.(fn(k)){l,i}{j};
                        end
                    end
                end
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Initialize field "e" (error): fields "e.a", "e.f" and "e.p".
        function [e] = Initialize_e(nc,ns,Nc,Nf)
            for i = ["a","f","p"]
                for j = 1:ns-1
                    %  > Error distribution.
                    if i ~= "f"
                        e.(i){j}.c.c_abs = zeros(Nc,1);
                        e.(i){j}.t.c     = zeros(Nc,1);
                        e.(i){j}.t.c_abs = zeros(Nc,1);
                        for k = ["x","y"]
                            e.(i){j}.t.f_abs.(k) = zeros(Nf,nc(1)+1);
                        end
                    else
                        e.(i){j}.c_abs = zeros(Nc,1);
                        for k = ["x","y"]
                            e.(i){j}.f_abs.(k) = zeros(Nf,nc(1)+1);
                        end
                    end
                    %  > Error norms.
                    e.(i){j}.n_abs.c = zeros(3,1);
                    if i ~= "f"
                        e.(i){j}.n_abs.t.c = zeros(3,1);
                        for k = ["x","y"]
                            e.(i){j}.n_abs.t.f.(k) = zeros(3,nc(1)+1);
                        end
                    else
                        for k = ["x","y"]
                            e.(i){j}.n_abs.f.(k) = zeros(3,nc(1)+1);
                        end
                    end
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update 'e.(...)' field (error).
        function [e] = Update_e(inp,msh,e,m,x)
            %  > Auxiliary variables.
            i      = 1;
            
            %  > #1: Update field 'e.a'.
            e.a{i} = B2_2D.Update_ea(msh,e.a{i},m{i},x{i});
        end
        %  > 1.2.1. -------------------------------------------------------
        %  > Update 'e.a(...)' field (analytic cell/face truncation error distribution/norms).
        function [e] = Update_ea(msh,e,m,x)
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
            fn = fieldnames(e.t.f_abs);
            l  = 1;
            for i = 1:msh.f.Nf
                c       = msh.f.ic  {i}(l);
                Sf(i,:) = msh.c.f.Sf{c}(msh.c.f.if(c,:) == i,:);
                for j = 1:n-1
                    for k = 1:n
                        e.t.f_abs.(fn{j})(i,k) = abs(e_t_fx{k}(i,j).*Sf(i,j));
                    end
                end
            end
            %  > Error norms.
            for j = 1:n-1
                e.n_abs.t.f.(fn{j}) = Tools_2.Set_n(e.t.f_abs.(fn{j}));
            end
            
            % >> c.
            %  > \tau_c.
            for i = 1:msh.c.Nc
                for j = 1:numel(msh.c.f.if(i,:))
                    e.t.c(i) = e.t.c(i)+msh.c.f.Sf{i}(j,:)*e_t_fx{n}(msh.c.f.if(i,j),:)';
                end
            end
            %  > \tau_c_abs.
            e.t.c_abs   = abs(e.t.c);
            %  > e_c_abs.
            %  > Equivalent to: ea.c.c_abs = f.av.c-x.nv.x.c.
            e.c.c_abs   = abs(m.At\e.t.c);
            %  > Error norms.
            e.n_abs.c   = Tools_2.Set_n(e.c.c_abs,msh.c.Volume);
            e.n_abs.t.c = Tools_2.Set_n(e.t.c_abs,msh.c.Volume);
        end
    end
end
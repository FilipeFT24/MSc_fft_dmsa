classdef B2_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "m" (matrices).
        function [m] = Initialize_3(nc,ns,Nc)
            for i = 1:ns
                for j = 1:nc(1)
                    for k = 1:nc(2)
                        m{i}.Ac{j}{k} = zeros(Nc);   %  > Ac{v}(x,y),Ac{g}(x,y).
                        m{i}.Bc{j}{k} = zeros(Nc,1); %  > Bc{v}(x,y),Bc{g}(x,y).
                    end
                end
                m{i}.At     = zeros(Nc);
                m{i}.Bt     = zeros(Nc,1);
                m{i}.nnz.Ac = zeros(1,nc(2));        %  > nnz(Ac(x),Ac(y)).
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
            %  > Rows (r) to be updated.
            for i = 1:numel(u.s)
                for j = 1:numel(u.s{i})
                    uc{i,j} = RunLength(sort([msh.f.ic{u.s{i}{j}}]));
                end
            end
            r         = RunLength(sort(cat(2,uc{:})));
            m.At(r,:) = 0;
            m.Bt(r,1) = f.st(r);
            
            %  > For each term (matrix)...
            for i = 1:numel(u.s)
                for j = 1:numel(u.s{i})
                    if ~isempty(u.s{i}{j})
                        %  > Reset...
                        m.Ac{i}{j}(uc{i,j},:) = 0;
                        m.Bc{i}{j}(uc{i,j},1) = 0;
                        %  > For each cell...
                        for k = uc{i,j}
                            %  > For each face...
                            for l = 1:numel(msh.c.f.if(k,:)), n = msh.c.f.if(k,l);
                                %  > Auxiliary variables.
                                o  = s.logical {n,i}{j};
                                a  = s.i       {n,i}{j}( o);
                                b  = s.i       {n,i}{j}(~o);
                                Sf = msh.c.f.Sf{k}     ( l,:);
                                %  > Ac (cell contributions).
                                m.Ac{i}{j}(k,a) = m.Ac{i}{j}(k,a)+Sf(j).*x.Tf_V{n,i}{j}(o);
                                %  > Bc (cell contributions).
                                if any(~o)
                                    i_bd            = arrayfun(@(x) find(f.bd.i == x),b);
                                    m.Bc{i}{j}(k,1) = m.Bc{i}{j}(k,1)-Sf(j)*x.Tf_V{n,i}{j}(~o)*f.bd.v(i_bd);
                                end
                            end
                        end
                    end
                    %  > Cumulative matrices.
                    m.At(r,:) = m.At(r,:)+m.Ac{i}{j}(r,:);
                    m.Bt(r,1) = m.Bt(r,1)+m.Bc{i}{j}(r,1);
                end
            end
            %  > nnz (x/y-contribution(s): v(x/y)+g(x/y)).
            for i = 1:numel(m.nnz.Ac)
                A       {i} = cellfun(@(x) x{:,i},m.Ac,'un',0);
                m.nnz.Ac(i) = nnz(sum(cat(3,A{i}{:}),3));
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
            flag = false;
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
                                x.vf.(fn(k)){l,i}{j}(~m,1) = f.bd.v(arrayfun(@(x) find(f.bd.i == x),b));
                            end
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
        function [e] = Initialize_e(nc,Nc,Nf)
            %  > Error distribution.
            for i = ["a","f","p"]
                if i ~= "f"
                    e.(i).c.c_abs = zeros(Nc,1);
                    e.(i).t.c     = zeros(Nc,1);
                    e.(i).t.c_abs = zeros(Nc,1);
                    e.(i).t.f_abs = zeros(Nf,nc(2)+1);
                else
                    e.(i) = cell(Nf,nc(2));
                end
            end
            %  > Error norms.
            for i = ["a","p"]
                e.(i).n_abs.c   = zeros(3,1);
                e.(i).n_abs.t.c = zeros(3,1);
                e.(i).n_abs.t.f = zeros(3,nc(2)+1);
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Update "e.(...)" field (error).
        function [e] = Update_e(inp,msh,e,m,s,x)
            %  > #1: Update field "e.a".
            e.a = B2_2D.Update_ea(msh,e.a,m,x);
        end
        %  > 2.2.1. -------------------------------------------------------
        %  > Update 'e.a(...)' field (analytic cell/face truncation error distribution/norms).
        function [e] = Update_ea(msh,e,m,x)
            %  > Reset...
            if any(e.t.c)
                e.t.c = zeros(msh.c.Nc,1);
            end
            
            % >> f.
            %  > \tau_f: {1}-\tau_f(\phi).
            %            {2}-\tau_f(\nabla\phi).
            %            {3}-\tau_f.
            n     = numel(x.xf.a);
            e_t_f = cell (1,n);
            for i = 1:n+1
                if i ~= n+1
                    e_t_f{i} = x.xf.a{i}-x.xf.x{i};
                else
                    e_t_f{i} = sum(cat(3,e_t_f{1:n}),3);
                end
            end
            %  > \tau_f_abs.
            k = 1;
            for i = 1:msh.f.Nf
                %  > Sf.
                c       = msh.f.ic  {i}(k);
                Sf(i,:) = msh.c.f.Sf{c}(msh.c.f.if(c,:) == i,:);
                %  > \tau_f{x,y,(t)} = [\tau_f*Sf(x),\tau_f*Sf(y),\tau_f*Sf(x,y)].
                for j = 1:n+1
                    if j ~= n+1
                        e.t.f_abs(i,j) = abs(e_t_f{n+1}(i,j)*Sf(i,j));
                    else
                        e.t.f_abs(i,j) = abs(e_t_f{n+1}(i,:)*Sf(i,:)');
                    end
                end
            end
            %  > Error norms.
            e.n_abs.t.f = src_Tools.Set_n(e.t.f_abs);
            
            % >> c.
            %  > \tau_c.
            for i = 1:numel(e.t.c)
                for j = 1:numel(msh.c.f.if(i,:))
                    e.t.c(i) = e.t.c(i)+msh.c.f.Sf{i}(j,:)*e_t_f{n}(msh.c.f.if(i,j),:)';
                end
            end
            %  > \tau_c_abs.
            e.t.c_abs   = abs(e.t.c);
            %  > e_c_abs.
            %  > Equivalent to: ea.c.c_abs = abs(f.av.c-x.nv.x.c).
            e.c.c_abs   = abs(m.At\e.t.c);
            %  > Error norms.
            e.n_abs.c   = src_Tools.Set_n(e.c.c_abs,msh.c.Volume);
            e.n_abs.t.c = src_Tools.Set_n(e.t.c_abs,msh.c.Volume);
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > 3.1.1. -------------------------------------------------------
        %  > Select faces for coarsening/refinement.
        function [u] = Update_u(inp,e,u)
            %  > Error treshold.
            trsh = inp.p_adapt.trsh.*max(e.t.f_abs(:,end));
            %  > Faces selected for refinement.
            fr_i = e.t.f_abs(:,1:end-1) > trsh;
            
            A = 2;
            for i = 1:size(u,2)
                u{i}.p{2}(fr_i) = u{i}.p{2}(fr_i)+A;
                u{i}.p;
                u{i}.s{1}{1} = [];
                u{i}.s{1}{2} = [];
                u{i}.s{2}{1} = find(fr_i(:,1));
                u{i}.s{2}{2} = find(fr_i(:,2));
            end
        end
        %  > 3.1.2. -------------------------------------------------------
        %  > Set stopping criterion/criteria.
        function [f] = Stop(inp,count,e)
            f = false;
            if count >= inp.p_adapt.nc, f = true;
                fprintf("Stopping criterion: max. number of cycles.\n");
            end
            if e(1)  <= inp.p_adapt.em, f = true;
                fprintf("Stopping criterion: min. error treshold (L1 norm).\n");
            end
        end
    end
end
classdef B2_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "m" (matrices).
        function [m] = Initialize_m(nc,ns,Nc)
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
        %  > Update matrix A and vector b.
        function [m] = Update_m(msh,f,m,s)
            %  > Rows (r) to be updated.
            for i = 1:numel(s.u.s)
                for j = 1:numel(s.u.s{i})
                    uc{i,j} = RunLength(sort([msh.f.ic{s.u.s{i}{j}}]));
                end
            end
            r         = RunLength(sort(cat(2,uc{:})));
            m.At(r,:) = 0;
            m.Bt(r,1) = f.st(r);
            
            %  > For each term...
            for i = 1:numel(s.u.s)
                for j = 1:numel(s.u.s{i})
                    if ~isempty(s.u.s{i}{j})
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
                                m.Ac{i}{j}(k,a) = m.Ac{i}{j}(k,a)+Sf(j).*s.tfV{n,i}{j}{1}(o);
                                %  > bc (face contributions).
                                if any(~o) 
                                    i_bd            = arrayfun(@(x) find(f.bd.i == x),b);
                                    m.Bc{i}{j}(k,1) = m.Bc{i}{j}(k,1)-Sf(j)*s.tfV{n,i}{j}{1}(~o)*f.bd.v(i_bd);
                                end
                                %  > Add kf (scalar) to the RHS...
                                if ~isempty(s.tfV{n,i}{j}{2})
                                    m.Bc{i}{j}(k,1) = m.Bc{i}{j}(k,1)-Sf(j).*s.tfV{n,i}{j}{2};
                                end
                            end
                        end
                    end
                    %  > At and bt(cumulative/total matrix/vector).
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
        % >> 1.3. ---------------------------------------------------------
        %  > Update nodal solution/face values (multiplied by V).
        function [s] = Update_sx(f,m,s)
            %  > Update nodal solution.
            s.x.nv.x = func.backlash(m.At,m.Bt);

            %  > Check if empty...
            fn   = convertCharsToStrings(fieldnames(s.x.vf));
            fn_n = numel(fn);
            flag = false;
            if all(cellfun(@isempty,s.x.vf.a),'all')
                flag = true;
            end
            %  > Update only necessary/all entries for field "a"/"x".
            for i = 1:numel(s.u.s)
                for j = 1:numel(s.u.s{i})
                    for k = 1:fn_n
                        switch fn(k)
                            case "a"
                                if ~flag
                                    r = s.u.s{i}{j}';
                                else
                                    r = 1:size(s.x.xfV.a{i},1);
                                end
                            case "x"
                                r = 1:size(s.x.xfV.a{i},1);
                        end
                        for l = r
                            %  > Auxiliary variables.
                            n = s.logical{l,i}{j};
                            a = s.i      {l,i}{j}( n);
                            b = s.i      {l,i}{j}(~n);
                            %  > Cell value(s).
                            s.x.vf.(fn(k)){l,i}{j}(n,1) = s.x.nv.(fn(k))(s.i{l,i}{j}(n));
                            %  > Face value(s).
                            if any(~n)
                                s.x.vf.(fn(k)){l,i}{j}(~n,1) = f.bd.v(arrayfun(@(x) find(f.bd.i == x),b));
                            end
                            %  > "xfV".
                            if ~isempty(s.tfV{l,i}{j}{2})
                                s.x.xfV.(fn(k)){i}(l,j) = s.tfV{l,i}{j}{1}*s.x.vf.(fn(k)){l,i}{j}+s.tfV{l,i}{j}{2}; %  > w/  constraint(s).
                            else
                                s.x.xfV.(fn(k)){i}(l,j) = s.tfV{l,i}{j}{1}*s.x.vf.(fn(k)){l,i}{j};                  %  > w/o constraint(s).
                            end
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
        %  > Update field "e.(...)" (error).
        function [e] = Update_e(inp,msh,e,m,s)
            %  > #1: Update field "e.a".
            e.a = B2_2D.Update_ea(msh,e.a,m,s);
        end
        %  > 2.2.1. -------------------------------------------------------
        %  > Update 'e.a(...)' field (analytic cell/face truncation error distribution/norms).
        function [e] = Update_ea(msh,e,m,s)
            %  > Reset...
            if any(e.t.c)
                e.t.c = zeros(msh.c.Nc,1);
            end
            
            % >> f.
            %  > \tau_f: {1}-\tau_f(\phi).
            %            {2}-\tau_f(\nabla\phi).
            %            {3}-\tau_f.
            n     = numel(s.x.xfV.a);
            e_t_f = cell (1,n);
            for i = 1:n+1
                if i ~= n+1
                    e_t_f{i} = s.x.xfV.a{i}-s.x.xfV.x{i};
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
            e.n_abs.t.f = func.Set_n(e.t.f_abs);
            
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
            e.c.c_abs   = abs(func.backlash(m.At,e.t.c));
            %  > Error norms.
            e.n_abs.c   = func.Set_n(e.c.c_abs,msh.c.Volume);
            e.n_abs.t.c = func.Set_n(e.t.c_abs,msh.c.Volume);
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > 3.1.1. -------------------------------------------------------
        %  > Check stopping criterion/criteria and continue/stop...
        function [f] = Stop(inp,cnt,e)
            f = false;
            if cnt  >= inp.p.n, f = true;
                fprintf("Stopping criterion: max. number of cycles.\n");
            end
            if e(1) <= inp.p.e, f = true;
                fprintf("Stopping criterion: min. error treshold (L1 norm).\n");
            end
        end
        %  > 3.1.2. -------------------------------------------------------
        %  > Select faces for coarsening/refinement.
        function [s_u] = Update_u(inp,msh,e,s_i,s_logical,s_u)
            %  > Auxiliary variables.
            %    NOTE: All terms are discretized in the same manner.
            A    = 2;
            k(1) = 1; %  > Convection/diffusion term.
            k(2) = 1; %  > x/y term.
            R    = 1;
            
            %  > Set structure "fr" (refinement).
            fr = B2_2D.Set_fr(inp,msh,e.t.f_abs(:,end),k,s_i,s_logical,s_u.p);
            %  > Set structure "u".
            %  > For each term... 
            %    NOTE: Do not change (convection and diffusion treated in a unified manner).
            for i = 1:numel(s_u.p)
                %  > For each direction... (x=y).
                if inp.p.iso
                    for j = 1:numel(s_u.p{i})
                        s_u.p{i}{j}(fr,:) = s_u.p{i}{j}(fr,:)+A;
                        s_u.s{i}{j}       = sort(fr);
                    end
                end
            end
        end
        %  > 3.1.2.1. -----------------------------------------------------
        function [lev] = lev_i(inp_p_iso,msh_c_f_if,up_k)
            %  > Isotropic.
            if inp_p_iso
                k(3) = 1;
                for j = 1:numel(msh_c_f_if)
                    lev(j) = ceil(up_k(msh_c_f_if(j),k(3))./2);
                end
            end
        end
        %  > 3.1.2.2. -----------------------------------------------------
        %  > ...for refinement.
        function [fr] = Set_fr(inp,msh,e,k,s_i,s_logical,s_up)
            %  > Set structure "sc".
            for i = 1:size(s_i,1)
                sc{i,1} = s_i{i,k(1)}{k(2)}(s_logical{i,k(1)}{k(2)});
            end

            %  > Set "lev" (cell order level).
            for i = 1:msh.c.Nc
                lev(i,:) = B2_2D.lev_i(inp.p.iso,msh.c.f.if(i,:),s_up{k(1)}{k(2)});
            end
            %  > Select faces based on error treshold.
            f      = 1:msh.f.Nf;
            trsh.v = inp.p.trsh(2).*max(e);
            trsh.i = f(e > trsh.v);
            %  > Check if the selected faces are eligible for refinement...
            cnt = 1;
            for i = trsh.i
                c = [msh.f.ic{i}]; n = numel(c); elig = zeros(1,n);
                for j = 1:n, c_j = c(j);
                    %  > Increment/update...
                    l        = msh.c.f.if(c_j,:) == i;
                    inc      = zeros(1,numel(l));
                    inc  (l) = 1;
                    lev_p    = lev(c_j,:)+inc;
                    lev_e    = e(msh.c.f.if(c_j,:))';
                    %  > Minimum/maximum p.
                    m        = min(lev_p); m_j = lev_p == m;
                    M        = max(lev_p); M_j = lev_p == M;
                    %  > Check...
                    if M-m <= 1
                        elig(j) = 1;
                    else
                        break;
                    end
                end
                if all(elig)
                    fr(cnt,1) = i; cnt = cnt+1;
                end
            end
        end
    end
end
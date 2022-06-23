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
            %  > Rows to be updated (r).
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
                A       {i} = cellfun(@(x) x{:,i},m.Ac,'un',false);
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
                            s.x.vf.(fn(k)){l,i}{j}      = zeros  (numel(n),1);
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
            %  > #2: Update field "e.f".
            e.f = B2_2D.Update_ef(msh,e.f,s);
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
        %  > 2.2.2. -------------------------------------------------------
        %  > Update 'e.f(...)' field (fitting error distribution).
        %  > e = sqrt(W)*[\phi_s-D*b], where b = (DTWD)^{-1}DTW*\phi_s.
        function [e] = Update_ef(msh,e,s)
            for i = 1:size(s.D,1)
                for j = 1:size(s.D,2)
                    for k = 1:numel(s.D{i,j})
                        e{i,j}{k} = abs(sqrt(s.D{i,j}{k}.W)*(s.x.vf.x{i,j}{k}-s.D{i,j}{k}.D*(s.D{i,j}{k}.DTWD\s.D{i,j}{k}.DTW*s.x.vf.x{i,j}{k})));
                    end
                end
            end 
        end
 
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > 3.1.1. -------------------------------------------------------
        %  > Select faces for coarsening/refinement.
        function [u,v,flag] = Update_u(inp,msh,e,u)
            %  > Auxiliary variables.
            %    NOTE: All terms are discretized in the same manner.
            A    = 2;
            k(1) = 1; %  > Convection/diffusion term.
            k(2) = 1; %  > x/y term.

            if inp.p.iso
                k(3) = 1;
                for i = 1:msh.f.Nf
                    lev.f(i,1) = ceil(u.p{k(1)}{k(2)}(i,k(3))./2);
                    F    {i,1} = B2_2D.find_F(msh,i);
                    lev.F{i,1} = ceil(u.p{k(1)}{k(2)}(F{i,1},k(3))./2);
                end
            end
            %  > Set structures "f" and "v".
            [f,v] = B2_2D.Set_fv(inp,msh,e.t.f_abs(:,end),F,lev);
            %  > Set structure "u".
            %  > For each term...
            %    NOTE: Do not change (convection and diffusion treated in a unified manner).
            for i = 1:numel(u.p)
                %  > For each direction... (x=y).
                if inp.p.iso
                    for j = 1:numel(u.p{i})
                        if ~isempty(f.c)
                            u.p{i}{j}(f.c,:) = u.p{i}{j}(f.c,:)-A; 
                        end
                        if ~isempty(f.r)
                            u.p{i}{j}(f.r,:) = u.p{i}{j}(f.r,:)+A;
                        end
                        u.s{i}{j} = sort([f.c;f.r]);
                    end
                end
            end
            if ~isempty(f.r) || ~isempty(f.c)
                flag = false;
            else
                flag = true;
            end
        end
        %  > 3.1.2. -------------------------------------------------------
        %  > Auxiliary function #1.
        %  > Face "f" vertex neighbours (layer #1).
        function [F] = find_F(msh,f)
            F = RunLength(sort(func.setdiff(cat(1,msh.v.if{msh.f.iv(f,:)})',f)));
        end
        %  > 3.1.3. -------------------------------------------------------
        %  > Auxiliary function #2.
        %  > Check whether face "f" is ilegible for coarsening/refinement...
        function [logical] = is_elig(ch,lev_f,lev_F)
            logical = false;
            switch ch
                case 'c', f_check = lev_f >= lev_F;
                case 'r', f_check = lev_f <= lev_F;
                otherwise
                    return;
            end
            if all(f_check)
                logical = true;
            end
        end
        %  > 3.1.4. -------------------------------------------------------
        %  > Auxiliary function #3.
        %  > Check what faces need to be coarsened/refined to coarsen/refine face "f"...
        function [add] = add_f(ch,f,F,lev)
            switch ch
                case 'c'
                    for i = 1:numel(f)
                        add(i).f = F{f(i)}(lev.F{f(i)} > lev.f(f(i)));
                    end
                case 'r'
                    for i = 1:numel(f)
                        add(i).f = F{f(i)}(lev.F{f(i)} < lev.f(f(i)));
                    end
                otherwise
                    return;
            end
            for i = 1:numel(f)
                for j = 1:numel(add(i).f)
                    add(i).logical(j) = B2_2D.is_elig(ch,lev.f(add(i).f(j)),lev.F{add(i).f(j)});
                end
            end
        end
        % >> 3.2. ---------------------------------------------------------
        %  > 3.2.1. -------------------------------------------------------
        function [f,v] = Set_fv(inp,msh,e,F,lev)
            % >> Set structure "list" and select faces for coarsening/refinement.
            %  > "i" and "v".
            [list.v,list.i] = sort(e,'descend');
            %  > "logical".
            for i = 1:numel(lev.f)
               list.logical.c(i,1) = B2_2D.is_elig('c',lev.f(i),lev.F{i}); 
               list.logical.r(i,1) = B2_2D.is_elig('r',lev.f(i),lev.F{i});
            end
            list.logical.c(:,2) = lev.f > 1;                    %  > Eligible for selection: p(f) > 1.
            list.logical.r(:,2) = lev.f < ceil(inp.p.p_max./2); %  > Eligible for selection: p(f) < p_max.
            list.logical.r      = list.logical.r(list.i,:);
            
            % >> Select faces and check whether they can be automatically...
            %  > ...refined.
            eval.r = B2_2D.f_select(inp,'r',list);
            if any(~eval.r.logical)
                %  > Identify valid/invalid faces.
                i_f.r = eval.r.i(~eval.r.logical);
                v_f.r = eval.r.i( eval.r.logical); f.r = v_f.r;
                %  > Loop through invalid faces...
                for i = 1:numel(i_f.r)
                    add_r = B2_2D.add_nb('r',i_f.r(i),F,lev);
                    f.r   = [f.r;add_r'];
                end
                f.r = RunLength(sort(f.r));
            else
                f.r = eval.r.i;
            end
            v.r = eval.r.i;
            
            %  > ...coarsened.
            %  > Select neighbours of ["f.r","v.r"]...
            tv_r                   = RunLength(sort(reshape(msh.f.iv(RunLength(sort([f.r;v.r])),:),1,[])));
            tn_r                   = RunLength(sort(cat(1,msh.v.if{tv_r})));
            list.logical.c(tn_r,2) = false;
            %  > Update field "F" of structure "lev" (w/ faces to be refined, i.e. "f.r")...
            for i = 1:numel(f.r)
                v_r{i,1} = RunLength   (sort(reshape(msh.f.iv(f.r(i),:),1,[])));
                n_r{i,1} = RunLength   (sort(cat(1,msh.v.if{v_r{i,1}})));
                n_r{i,1} = func.setdiff(n_r{i,1},f.r(i));
                for j = 1:numel(n_r{i,1})
                    lev.F{n_r{i,1}(j)}(F{n_r{i,1}(j)} == f.r(i)) = ...
                        lev.F{n_r{i,1}(j)}(F{n_r{i,1}(j)} == f.r(i))+1;
                end
            end
            for i = func.setdiff(1:numel(lev.F),tn_r)
                for j = 1:numel(F{i})
                    if range(lev.F{F{i}(j)}) > 1
                        list.logical.c(i,2) = false;
                        break;
                    end
                end
            end
            list.logical.c = list.logical.c(list.i,:);

            eval.c = B2_2D.f_select(inp,'c',list);
            if any(~eval.c.logical)
                %  > Identify valid/invalid faces.
                i_f.c = eval.c.i(~eval.c.logical);
                v_f.c = eval.c.i( eval.c.logical); f.c = v_f.c;
                %  > Loop through invalid faces...
                for i = 1:numel(i_f.c)
                    add_c = B2_2D.add_nb('c',i_f.c(i),F,lev);
                    f.c   = [f.c;add_c'];
                end
                f.c = RunLength(sort(f.c));
            else
                f.c = eval.c.i;
            end
            v.c = eval.c.i;
        end
        %  > 3.2.2. -------------------------------------------------------
        function [eval] = f_select(inp,ch,list)
            switch ch
                case 'c'
                    %  > ...for coarsening.
                    n            = numel (list.v);
                    v            = list.v(n-ceil(inp.p.trsh(1)*n));
                    cond         = list.logical.c(:,2) & list.v <= v;
                    eval.i       = list.i        (cond);
                    eval.logical = list.logical.c(cond,1);
                case 'r'
                    %  > ...for refinement.
                    n            = numel (list.v);
                    v            = list.v(n-ceil(inp.p.trsh(2)*n));
                    cond         = list.logical.r(:,2) & list.v >= v;
                    eval.i       = list.i        (cond);
                    eval.logical = list.logical.r(cond,1);
                otherwise
                    return;
            end
        end
        %  > 3.2.3. -------------------------------------------------------
        function [add_all] = add_nb(ch,f,F,lev)
            %  > Add...
            f_add = B2_2D.add_f(ch,f,F,lev);
            %  > Assign to "add_all".
            if all(f_add.logical)
                add_all = f_add.f;
            else
                g = f_add.f(~f_add.logical);
                for i = 1:numel(g), j = g(i); fr{i} = [];
                    while 1
                        %  > Add...
                        add_g = B2_2D.add_f(ch,j,F,lev);
                        %  > Check...
                        if all([add_g(:).logical])
                            fr{i} = [fr{i},[add_g(:).f]];
                            break;
                        else
                            h     =    [add_g(:).f];
                            fr{i} = h( [add_g(:).logical]);
                            j     = h(~[add_g(:).logical]);
                        end
                    end
                end
                add_all = [f_add.f(f_add.logical),[fr{:}]];
            end
        end
        % >> 3.3. ---------------------------------------------------------
        %  > Check stopping criteria.
        function [stop] = Stop(inp,flag,ec,tau_f)
            %  > Auxiliary variables.
            n    = numel(tau_f);
            stop = false(1,3);
            
            %  > (1).
            if n     >  inp.p.N, stop(1) = true;
                fprintf("Stopping criterion: max. number of cycles.\n");
            end
            %  > (2).
            if ec(1) <= inp.p.e, stop(2) = true;
                fprintf("Stopping criterion: min. error treshold (L1 norm).\n");
            end
            %  > (3).
            if n > inp.p.n
                if all(tau_f(n-inp.p.n:n-1) < tau_f(n-inp.p.n+1:n)), stop(3) = true;
                    fprintf("Stopping criterion: increasing tau_f (L1 norm). \n%d previous iterations have been removed).",inp.p.n);
                end
            end
            %  > (4)
            if flag, stop(4) = true;
                fprintf("Stopping criterion: couldn't coarsen/refine any further.\n");
            end
        end
    end
end
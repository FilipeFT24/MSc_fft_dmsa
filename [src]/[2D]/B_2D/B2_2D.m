classdef B2_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "m" (matrices).
        function [m] = Initialize_m(ns,Nc)
            for i = 1:ns
                m{i}.A     = zeros(Nc);
                m{i}.b     = zeros(Nc,1);
                m{i}.nnz.A = 0;
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update matrix A and vector b.
        function [m] = Update_m(msh,f,m,s)
            %  > Rows to be updated (r).
            r        = RunLength(sort([msh.f.ic{s.u.s}]));
            m.A(r,:) = 0;
            m.b(r,1) = f.st(r);
            
            %  > For each cell...
            for i = r
                %  > For each face...
                for j  = 1:numel(msh.c.f.if(i,:)), k = msh.c.f.if(i,j);
                    %  > Auxiliary variables.
                    o  = s.logical {k};
                    a  = s.i       {k}( o);
                    b  = s.i       {k}(~o);
                    Sf = msh.c.f.Sf{i}( j,:);
                    %  > For each term...
                    for l = 1:size(s.x.s.tfV,2)
                        %  > A.
                        m.A(i,a) = m.A(i,a)+Sf*s.x.s.tfV{k,l}{1}(:,o);
                        %  > b.
                        if any(~o)
                            m.b(i,1) = m.b(i,1)-Sf*s.x.s.tfV{k,l}{1}(:,~o)*f.bd.v(func.find_c(f.bd.i,b));
                        end
                        %  > Add kf (scalar) to the RHS...
                        if ~isempty(s.x.s.tfV{k,l}{2})
                            m.b(i,1) = m.b(i,1)-Sf*s.x.s.tfV{k,l}{2};
                        end
                    end
                end
            end
            m.nnz.A = nnz(m.A);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Update nodal solution/face values (multiplied by V).
        function [s] = Update_sx(f,m,s)
            %  > Auxiliary variables.
            fn = convertCharsToStrings(fieldnames(s.x.x.vf));
            if any(~cellfun(@isempty,s.x.x.vf.a))
                c = s.u.s';
            else
                c = 1:size(s.x.x.vf.a,1);
            end
            
            %  > Update nodal solution.
            s.x.x.nv.x = m.A\m.b;
            %  > Update face coefficients/face values.
            for i = 1:numel(fn)
                for j = c
                    % >> Update nodal values used to fit face "j"...
                    %  > Auxiliary variables.
                    n = s.logical{j};
                    a = s.i      {j}( n);
                    b = s.i      {j}(~n);
                    %  > Cell value(s).
                    s.x.x.vf.(fn(i)){j}      = zeros(numel(n),1);
                    s.x.x.vf.(fn(i)){j}(n,1) = s.x.x.nv.(fn(i))(s.i{j}(n));
                    %  > Face value(s).
                    if any(~n)
                        s.x.x.vf.(fn(i)){j}(~n,1) = f.bd.v(func.find_c(f.bd.i,b));
                    end
                    % >> Update coefficients/face values...
                    for k = 1:size(s.x.s.tfV,2)
                        if ~isempty(s.x.s.tfV{j,k}{2})
                            s.x.x.xfV.(fn(i)){k}(j,:) = s.x.s.tfV{j,k}{1}*s.x.x.vf.(fn(i)){j}+s.x.s.tfV{j,k}{2}; %  > w/  constraint(s).
                        else
                            s.x.x.xfV.(fn(i)){k}(j,:) = s.x.s.tfV{j,k}{1}*s.x.x.vf.(fn(i)){j};                   %  > w/o constraint(s).
                        end
                    end
                end
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Initialize field "e" (error): fields "e.a", "e.f" and "e.p".
        function [e] = Initialize_e(Nc,Nf)
            %  > Error distribution.
            for i = ["a","f","p"]
                if i ~= "f"
                    e.(i).c.c_abs = zeros(Nc,1);
                    e.(i).t.c     = zeros(Nc,1);
                    e.(i).t.c_abs = zeros(Nc,1);
                    e.(i).t.f_abs = zeros(Nf,1);
                else
                    e.(i) = cell(Nf,1);
                end
            end
            %  > Error norms.
            for i = ["a","p"]
                e.(i).n_abs.c   = zeros(3,1);
                e.(i).n_abs.t.c = zeros(3,1);
                e.(i).n_abs.t.f = zeros(3,1);
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
        %  > Update "e.a(...)" field (analytic cell/face truncation error distribution/norms).
        function [e] = Update_ea(msh,e,m,s)
            %  > Reset...
            if any(e.t.c)
                e.t.c = zeros(msh.c.Nc,1);
            end
            
            % >> f.
            %  > \tau_f: {1}-\tau_f(\phi).
            %            {2}-\tau_f(\nabla\phi).
            %            {3}-\tau_f.
            n     = numel(s.x.x.xfV.a);
            e_t_f = cell (1,n);
            for i = 1:n+1
                if i ~= n+1
                    e_t_f{i} = s.x.x.xfV.a{i}-s.x.x.xfV.x{i};
                else
                    e_t_f{i} = sum(cat(3,e_t_f{1:n}),3);
                end
            end
            %  > \tau_f_abs.
            k = 1;
            for i = 1:msh.f.Nf
                c              = msh.f.ic  {i}(k);
                Sf       (i,:) = msh.c.f.Sf{c}(msh.c.f.if(c,:) == i,:);
                e.t.f_abs(i,1) = abs(e_t_f{n+1}(i,:)*Sf(i,:)');
            end
            %  > Error norms.
            e.n_abs.t.f = func.Set_n(e.t.f_abs);
            
            % >> c.
            %  > \tau_c.
            for i = 1:msh.c.Nc
                e.t.c(i,1) = sum(sum(msh.c.f.Sf{i}'.*e_t_f{n+1}(msh.c.f.if(i,:),:)',1),2);
            end
            %  > \tau_c_abs.
            e.t.c_abs   = abs(e.t.c);
            %  > e_c_abs.
            e.c.c_abs   = abs(s.x.x.nv.a-s.x.x.nv.x);
            %           = abs(m.A\e.t.c);
            %  > Error norms.
            e.n_abs.c   = func.Set_n(e.c.c_abs,msh.c.Volume);
            e.n_abs.t.c = func.Set_n(e.t.c_abs,msh.c.Volume);
        end
        %  > 2.2.2. -------------------------------------------------------
        %  > Update "e.f(...)" field (fitting error distribution).
        %  > e = sqrt(W)*[\phi_s-D*b], where b = (DTWD\DTW)*\phi_s.
        function [e] = Update_ef(msh,e,s)
            for i = 1:size(s.D,1)
                b{i} = s.D{i}.DTWD\s.D{i}.DTW*s.x.x.vf.x{i};
                e{i} = abs(sqrt(s.D{i}.W)*(s.x.x.vf.x{i}-s.D{i}.D*b{i}));
            end
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > 3.1.1. -------------------------------------------------------
        %  > Select faces for coarsening/refinement.
        function [u,v,flag] = Update_u(inp,msh,e,u)
            %  > Set structures "f" and "v".
            if inp.p.iso
                j = 1;
                for i = 1:size(u.p,1)
                    lev.f(i,1) = ceil(u.p(i,j)./2);
                    F    {i,1} = B2_2D.find_F(msh,i);
                    lev.F{i,1} = ceil(u.p(F{i,1},j)./2);
                end
            end
            [f,v] = B2_2D.Set_fv(inp,msh,e.t.f_abs,F,lev);
            %  > Set structure "u".
            if isempty(f.r) && isempty(f.c)
                flag       = true;
            else
                A          = 2;
                u.p(f.c,:) = u.p(f.c,:)-A;
                u.p(f.r,:) = u.p(f.r,:)+A;
                u.s        = sort([f.c;f.r]);
                flag       = false;
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
                    cond         = list.logical.c(:,2) & list.v <= round(v,10);
                    eval.i       = list.i        (cond);
                    eval.logical = list.logical.c(cond,1);
                case 'r'
                    %  > ...for refinement.
                    n            = numel (list.v);
                    v            = list.v(n-ceil(inp.p.trsh(2)*n));
                    cond         = list.logical.r(:,2) & list.v >= round(v,10);
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
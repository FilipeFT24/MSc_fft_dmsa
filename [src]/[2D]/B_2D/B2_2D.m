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
        %  > Check stopping criteria.
        %  > 3.1.1. -------------------------------------------------------
        %  > Auxiliary function #1.
        function [f] = Stop_1(inp,i,e)
            f = false;
            if i    >  inp.p.N, f = true;
                fprintf("Stopping criterion: max. number of cycles.\n");
            end
            if e(1) <= inp.p.e, f = true;
                fprintf("Stopping criterion: min. error treshold (L1 norm).\n");
            end
        end
        %  > 3.1.2. -------------------------------------------------------
        %  > Auxiliary function #2.
        function [f] = Stop_2(e,n)
            f = false;
            if all(e(end-n:end-1) < e(end-n+1:end))
                f = true;
                fprintf("Stopping criterion: increasing (face) truncation error (L_infinity norm).\n");
            end
        end
        % >> 3.2. ---------------------------------------------------------
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
        %  > 3.2.1. -------------------------------------------------------
        %  > 3.2.1.1. -----------------------------------------------------
        %  > ...for refinement.
        function [fr] = Set_fr(inp,msh,e,k,s_i,s_logical,s_up)
            %  > Auxiliary variables.
            Nc           = msh.c.Nc;
            Nf           = msh.f.Nf;
            list.logical = true(Nf,2);
            if inp.p.iso
                k(3) = 1;
                for i = 1:Nf
                    lev_f(i,1) = ceil(s_up{k(1)}{k(2)}(i,k(3))./2);
                end
            end
            
            % >> 1. Sort by "\tau_f".
            [list.v,list.i] = sort(e,'descend');
            %  > Check...
            for i = 1:Nf
               list.logical(i) = B2_2D.ref_check_1(inp,msh,i,lev_f);
            end
            list.logical = list.logical(list.i);
            % >> 2. Select faces...
            cond   = list.v > inp.p.trsh(2).*list.v(1);
            elig_f = list.logical(cond);
            elig_i = list.i      (cond);
            %  > Check...
            if any(~elig_f)
                %  > Identify eligible/ilegible faces.
                fr = elig_i( elig_f);
                gr = elig_i(~elig_f);
                %  > Loop through ilegible faces...
                for i = 1:numel(gr)
                    fr_add = B2_2D.ref_check_2(inp,msh,gr(i),lev_f);
                    fr     = [fr;fr_add'];
                end
                fr = RunLength(sort(fr));
            else
                fr = list.i(elig_f);
            end
        end
        %  > 3.2.1.2. -----------------------------------------------------
        %  > ...for refinement (auxiliary function #1).
        function [logical] = ref_check_1(inp,msh,f,lev_f)
            %  > Faces that contain face "f"'s vertices.
            fv = B2_2D.find_fv(msh,f);
            %  > Check whether face "f" is ilegible for refinement...
            if lev_f(f)+1 <= ceil(inp.p.p_max./2)
                f_check = lev_f(f) <= lev_f(fv);
                if all(f_check)
                    logical = true;
                else
                    logical = false;
                end
            else
                logical = false;
            end
        end
        %  > 3.2.1.3. -----------------------------------------------------
        %  > ...for refinement (auxiliary function #2).
        function [fr_all] = ref_check_2(inp,msh,f,lev_f)
            %  > Add...
            add_f = B2_2D.Add_f(inp,msh,f,lev_f);
            %  > Assign to "fr".
            if all(add_f.logical)
                fr_all = add_f.f;
            else
                g = add_f.f(~add_f.logical);
                for i = 1:numel(g), h = g(i); fr{i} = [];
                    while 1
                        %  > Add...
                        add_g = B2_2D.Add_f(inp,msh,h,lev_f);
                        %  > Check...
                        if all([add_g(:).logical])
                            fr{i} = [fr{i},[add_g(:).f]];
                            break;
                        else
                            fr{i} = add_g.f( [add_g(:).logical]);
                            h     = add_g.f(~[add_g(:).logical]);
                        end
                    end
                end
                fr_all = [add_f.f(add_f.logical),[fr{:}]];
            end
        end
        %  > 3.2.1.4. -----------------------------------------------------
        %  > ...for refinement (auxiliary function #3).
        function [fv] = find_fv(msh,f)
            fv = RunLength(sort(func.setdiff(cat(1,msh.v.if{msh.f.iv(f,:)})',f)));
        end
        %  > 3.2.1.5. -----------------------------------------------------
        %  > ...for refinement (auxiliary function #4).
        function [add] = Add_f(inp,msh,f,lev_f)
            for i = 1:numel(f)
                %  > Faces that contain face "f"'s vertices.
                fv{i} = B2_2D.find_fv(msh,f(i));
                %  > Check...
                add(i).f = fv{i}(lev_f(fv{i}) < lev_f(f(i)));
                for j = 1:numel(add(i).f)
                    add(i).logical(j) = B2_2D.ref_check_1(inp,msh,add(i).f(j),lev_f);
                end
            end
        end
    end
end
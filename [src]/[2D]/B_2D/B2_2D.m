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
                m.At(r,:) = m.At(r,:)+m.Ac{i}(r,:);
                m.Bt(r,1) = m.Bt(r,1)+m.Bc{i}(r,1);
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
        function [e] = Initialize_e(nc,Nc,Nf)
            %  > Error distribution.
            for i = ["a","f","p"]
                if i ~= "f"
                    e.(i).c.c_abs    = zeros(Nc,1);
                    e.(i).t.c        = zeros(Nc,1);
                    e.(i).t.c_abs    = zeros(Nc,1);
                    for j = 1:numel(nc)+1
                        if j ~= numel(nc)+1
                            e.(i).t.f_abs{j} = zeros(Nf,nc(j)+1);   %  > v,g,(t) and x,y,(t).
                        else
                            e.(i).t.f_abs{j} = zeros(Nf,sum(nc)+1); %  > vx,vy,gx,gy,(t).
                        end
                    end
                else
                    e.(i).c_abs = zeros(Nc,1);
                    for j = 1:nc(1)
                        e.(i).f_abs{j} = zeros(Nf,nc(2));   %  > v(x,y),g(x,y).
                    end
                end
            end
            %  > Error norms.
            for i = ["a","f","p"]
                e.(i).n_abs.c = zeros(3,1);
                if i ~= "f"
                    e.(i).n_abs.t.c    = zeros(3,1);
                    for j = 1:numel(nc)+1
                        if j ~= numel(nc)+1
                            e.(i).n_abs.t.f{j} = zeros(3,nc(j)+1);   %  > v,g,(t) and x,y,(t).
                        else
                            e.(i).n_abs.t.f{j} = zeros(3,sum(nc)+1); %  > vx,vy,gx,gy,(t).
                        end
                    end
                else
                    for j = 1:nc(1)
                        e.(i).n_abs.f{j} = zeros(3,nc(2)); %  > v(x,y),g(x,y).
                    end
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update 'e.(...)' field (error).
        function [e] = Update_e(inp,msh,e,m,x)
            %  > #1: Update field 'e.a'.
            i   = 1;
            e.a = B2_2D.Update_ea(msh,e.a,m{i},x{i});
        end
        %  > 1.2.1. -------------------------------------------------------
        %  > Update 'e.a(...)' field (analytic cell/face truncation error distribution/norms).
        function [e] = Update_ea(msh,e,m,x)
            %  > Initialize...
            e.t.c = zeros(msh.c.Nc,1);
            
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
            e_t_f_3 = cat(2,e_t_f{1:n});
            %  > \tau_f_abs.
            l = 1;
            for i = 1:msh.f.Nf
                %  > Sf.
                c       = msh.f.ic  {i}(l);
                Sf(i,:) = msh.c.f.Sf{c}(msh.c.f.if(c,:) == i,:);
                %  > #1: v,g,(t).
                %  > [\tau_f\phi(x,y),\tau_f\nabla\phi(x,y),\tau_f(x,y)]*Sf(x,y).
                for j = 1:size(e_t_f,2)
                    e.t.f_abs{1}(i,j) = abs(e_t_f{j}(i,:)*Sf(i,:)');
                end
                %  > #2: x,y,(t).
                %  > [\tau_f*Sf(x),\tau_f*Sf(y),\tau_f*Sf(x,y)].
                for j = 1:size(e_t_f{end},2)+1
                    if j ~= size(e_t_f{end},2)+1
                        e.t.f_abs{2}(i,j) = abs(e_t_f{end}(i,j)*Sf(i,j));
                    else
                        e.t.f_abs{2}(i,j) = abs(e_t_f{end}(i,:)*Sf(i,:)');
                    end
                end
                %  > #3: vx,vy,gx,gy,(t).
                %  > [\tau_f\phi(x)*Sf(x),\tau_f\phi(y)*Sf(y),\tau_f\nabla\phi(x)*Sf(x),\tau_f\nabla\phi(y)*Sf(y),\tau_f(x,y)*Sf(x,y)].
                %  > Multiply by Sf(x)...
                e.t.f_abs{3}(i,[1,3]) = e_t_f_3(i,[1,3]).*Sf(i,1);
                %  > Multiply by Sf(y)...
                e.t.f_abs{3}(i,[2,4]) = e_t_f_3(i,[2,4]).*Sf(i,2);
            end
            %  > Sum contributions (3)...
            e.t.f_abs{3}(:,end) = sum(e.t.f_abs{3}(:,1:end-1),2);
            e.t.f_abs{3}        = abs(e.t.f_abs{3});
            %  > Error norms.
            for i = 1:numel(e.n_abs.t.f)
                e.n_abs.t.f{i} = Tools_2.Set_n(e.t.f_abs{i});
            end
            
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
            e.n_abs.c   = Tools_2.Set_n(e.c.c_abs,msh.c.Volume);
            e.n_abs.t.c = Tools_2.Set_n(e.t.c_abs,msh.c.Volume);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > 1.3.1. -------------------------------------------------------
        %  > Select faces for coarsening/refinement.
        function [u] = Update_u(inp,e,u)
            A = 2; %  > Method's order.
            k = 2; %  > Choice: #1: v,g,(t).
                   %            #2: x,y,(t).
                   %            #3: vx,vy,gx,gy,(t).
            
            %  > Error treshold.
            trsh = inp.p_adapt.trsh.*max(e.t.f_abs{k}(:,end));
            %  > Faces selected for refinement.
            fr_i = e.t.f_abs{k}(:,1:end-1) > trsh;
            
            for i = 1:size(u,2)
                u{i}.p{2}(fr_i) = u{i}.p{2}(fr_i)+A;
                u{i}.p;
                
                u{i}.s{1}{1} = [];
                u{i}.s{1}{2} = [];
                u{i}.s{2}{1} = find(fr_i(:,1));
                u{i}.s{2}{2} = find(fr_i(:,2));
            end
        end
        %  > 1.3.2. -------------------------------------------------------
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
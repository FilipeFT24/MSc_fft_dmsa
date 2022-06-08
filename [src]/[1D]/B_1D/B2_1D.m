classdef B2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "m" (matrices).
        function [m] = Initialize_m(nc,ns,Nc)
            for i = 1:ns
                for j = 1:nc
                    m{i}.Ac{j} = zeros(Nc);   %  > Ac{v},Ac{g}.
                    m{i}.bc{j} = zeros(Nc,1); %  > bc{v},bc{g}.
                end
                m{i}.At     = zeros(Nc);
                m{i}.bt     = zeros(Nc,1);
                m{i}.nnz.Ac = zeros(1,nc);    %  > nnz(Ac).
                m{i}.nnz.At = 0;
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update matrix A and vector b.
        function [m] = Update_m(msh,f,m,s)
            %  > Rows(r) to be updated.
            for i = 1:numel(s.u.s)
                uc{i} = RunLength([msh.f.ic{s.u.s{i}}]);
            end
            r         = RunLength(cat(2,uc{:}));
            m.At(r,:) = 0;
            m.bt(r,1) = f.st(r);
            
            %  > For each term...
            for i = 1:numel(s.u.s)
                %  > Reset...
                m.Ac{i}(uc{i},:) = 0;
                m.bc{i}(uc{i},1) = 0;
                %  > For each cell...
                for j = uc{i}
                    %  > For each face...
                    for k = 1:numel(msh.c.f.if(j,:)), n = msh.c.f.if(j,k);
                        %  > Auxiliary variables.
                        l = s.logical{n,i};
                        a = s.i      {n,i}( l);
                        b = s.i      {n,i}(~l);
                        %  > Ac(cell contributions).
                        m.Ac{i}(j,a) = m.Ac{i}(j,a)+msh.c.f.Sf(j,k)*s.tfV{n,i}(l);
                        %  > bc(cell contributions).
                        if any(~l)
                            i_bd         = arrayfun(@(x) find(f.bd.i == x),b);
                            m.bc{i}(j,1) = m.bc{i}(j)-msh.c.f.Sf(j,k)*s.tfV{n,i}(~l)*f.bd.v(i_bd);
                        end
                    end
                end
                %  > At and bt(cumulative/total matrix/vector).
                m.At = m.At+m.Ac{i};
                m.bt = m.bt+m.bc{i};
            end
            %  > nnz.
            for i = 1:numel(s.u.s)
                m.nnz.Ac(i) = nnz(m.Ac{i});
            end
            m.nnz.At = nnz(m.At);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Update nodal solution/face values (multiplied by V).
        function [s] = Update_sx(f,m,s)
            %  > Update nodal solution.
            s.x.nv.x = m.At\m.bt;
            
            %  > Check if empty...
            fn   = convertCharsToStrings(fieldnames(s.x.vf));
            fn_n = numel(fn);
            flag = false;
            if all(cellfun(@isempty,s.x.vf.a),'all')
                flag = true;
            end
            %  > Update only necessary/all entries for fields "a"/"x".
            for i = 1:numel(s.u.s)
                for j = 1:fn_n
                    switch fn(j)
                        case "a"
                            if ~flag
                                r = s.u.s{i}';
                            else
                                r = 1:size(s.x.xfV.a,1);
                            end
                        case "x"
                            r = 1:size(s.x.xfV.a,1);
                    end
                    for k = r
                        %  > Auxiliary variables.
                        m = s.logical{k,i};
                        a = s.i      {k,i}( m);
                        b = s.i      {k,i}(~m);
                        %  > Cell value(s).
                        s.x.vf.(fn(j)){k,i}(m,1) = s.x.nv.(fn(j))(s.i{k,i}(m));
                        %  > Face value(s).
                        if any(~m)
                            s.x.vf.(fn(j)){k,i}(~m,1) = f.bd.v(arrayfun(@(x) find(f.bd.i == x),b));
                        end
                        %  > "xfV".
                        s.x.xfV.(fn(j))(k,i) = s.tfV{k,i}*s.x.vf.(fn(j)){k,i};
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
                    e.(i).t.f_abs = zeros(Nf,nc+1);
                else
                    e.(i) = cell(Nf,nc);
                end
            end
            %  > Error norms.
            for i = ["a","p"]
                e.(i).n_abs.c   = zeros(3,1);
                e.(i).n_abs.t.c = zeros(3,1);
                e.(i).n_abs.t.f = zeros(3,nc+1);
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Update field "e.(...)" (error).
        function [e] = Update_e(inp,msh,e,m,s)
            %  > #1: Update field "e.a".
            e.a = B2_1D.Update_ea(msh,e.a,m,s);
        end
        %  > 2.2.1. -------------------------------------------------------
        %  > Update field "e.a(...)" (analytic cell/face truncation error distribution/norms).
        function [e] = Update_ea(msh,e,m,s)
            %  > Reset...
            if any(e.t.c)
                e.t.c = zeros(msh.c.Nc,1);
            end
            
            % >> f.
            %  > \tau_f: {1}-\tau_f(\phi).
            %            {2}-\tau_f(\nabla\phi).
            %            {3}-\tau_f.
            n = size(s.x.xfV.a,2);
            for i = 1:n+1
                if i ~= n+1
                    e_t_f(:,i) = s.x.xfV.a(:,i)-s.x.xfV.x(:,i);
                else
                    e_t_f(:,i) = sum(e_t_f(:,1:n),2);
                end
            end
            %  > \tau_f_abs.
            e.t.f_abs = abs(e_t_f);
            %  > Error norms.
            e.n_abs.t.f = src_Tools.Set_n(e.t.f_abs);
            
            % >> c.
            %  > \tau_c.
            for i = 1:numel(e.t.c)
                for j = 1:numel(msh.c.f.if(i,:))
                    e.t.c(i) = e.t.c(i)+msh.c.f.Sf(i,j)*e_t_f(msh.c.f.if(i,j),n+1)';
                end
            end
            %  > \tau_c_abs.
            e.t.c_abs   = abs(e.t.c);
            %  > e_c_abs.
            e.c.c_abs   = abs(m.At\e.t.c);
            %  > Error norms.
            e.n_abs.c   = src_Tools.Set_n(e.c.c_abs,msh.c.Volume);
            e.n_abs.t.c = src_Tools.Set_n(e.t.c_abs,msh.c.Volume);
        end
        %  > 2.2.2. -------------------------------------------------------
        %  > Update field "e.p(...)" (predicted cell/face truncation error distribution/norms).
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > 3.1.1. -------------------------------------------------------
        %  > Select faces for coarsening/refinement.
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
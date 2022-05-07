classdef B2_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "m" (matrices).
        function [m] = Initialize_3(nc,ns,Nc)
            for i = 1:ns
                for j = 1:nc
                    m{i}.Ac{j} = zeros(Nc);
                    m{i}.Bc{j} = zeros(Nc,1);
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
            for i = 1:size(u.s,2)
                cf{i} = [msh.f.ic{u.s{i}}]';
            end
            r         = unique(cat(1,cf{:}));
            m.At(r,:) = 0;
            m.Bt(r)   = f.st(r);
   
            %  > Loop through direction(s)...
            for i = 1:size(u.s,2)
                %  > For each cell...
                for j = r'
                    %  > Loop through cell "j"'s faces.
                    for k = 1:numel(msh.c.f.if(j,:))
                        %  > Face index (f_jk).
                        f_jk         = msh.c.f.if(j,k);
                        %  > Cell/face indices used to fit "f_jk".
                        l            = s.logical{f_jk,i};
                        a            = s.i      {f_jk,i}( l);
                        b            = s.i      {f_jk,i}(~l); b = sort(b);
                        %  > Ac (cell contributions).
                        m.Ac{i}(j,a) = m.Ac{i}(j,a)+msh.c.f.Sf{j}(k,:)*x.Tf_V{f_jk,i}(:,l);
                        %  > Bc (cell contributions).
                        if any(~l)
                            bd_v       = f.bd.v (ismembc(f.bd.i,b));
                            m.Bc{i}(j) = m.Bc{i}(j)-msh.c.f.Sf{j}(k,:)*x.Tf_V{f_jk,i}(:,~l)*bd_v;
                        end
                    end
                end
                %  > At and Bt (cumulative/total matrices).
                m.At = m.At+m.Ac{i};
                m.Bt = m.Bt+m.Bc{i};
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
                    e.(i){j}.t.c       = zeros(Nc,1);
                    e.(i){j}.t.c_abs   = zeros(Nc,1);
                    e.(i){j}.t.f       = cell (1 ,nc+1);
                    e.(i){j}.t.f_abs   = cell (1 ,nc+1);
                    for k = 1:size(e.(i){j}.t.f,2)
                        e.(i){j}.t.f    {k} = zeros(Nf,nc);
                        e.(i){j}.t.f_abs{k} = zeros(Nf,nc);
                    end
                    e.(i){j}.t.n.c     = zeros(1 ,3);
                    e.(i){j}.t.n_abs.c = zeros(1 ,3);
                    e.(i){j}.t.n.f     = zeros(3 ,nc+1);
                    e.(i){j}.t.n_abs.f = zeros(3 ,nc+1);
                    e.(i){j}.c.c       = zeros(Nc,1);
                    e.(i){j}.c.c_abs   = zeros(Nc,1);
                    e.(i){j}.c.n       = zeros(1 ,3);
                    e.(i){j}.c.n_abs   = zeros(1 ,3);
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update error.
        function [e] = Update_e(inp,msh,e,f,m,s,u,x)
            
            j = 1;
            e.a{j} = B2_2D.Set_1_e(msh,e.a{j},inp.c,x{j}.xf.a,x{j}.nv.a.f);
            
            e.a{j} = B2_2D.Set_2_e(e.a{j},m{j},msh.c.Volume);
            
        end
        % >> 1.3. ---------------------------------------------------------
        %  > 1.3.1. -------------------------------------------------------
        %  > Auxiliary function #1.
        %  > Compute remaining error fields (based on the convective/diffusive facial components).
        function [e] = Set_1_e(msh,e,v,x,y)
            %  > \tau_f(\phi)_(x,y) and \tau(\nabla\phi)_(x,y).
            [~,a] = size(e.t.f);
            for i = 1:a-1
                e.t.f{i} = v(i,:).*(y{i}-x{i});
            end
            %  > \tau_f.
            e.t.f{a} = sum(cat(3,e.t.f{1:a-1}),3);
            %  > \tau_c.
            for i = 1:size(e.t.c,1)
                for j = 1:numel(msh.c.f.if(i,:))
                    e.t.c(i) = e.t.c(i)+msh.c.f.Sf{i}(j,:)*e.t.f{a}(msh.c.f.if(i,j),:)';
                end
            end
        end
        %  > 1.3.2. -------------------------------------------------------
        %  > Auxiliary function #2.
        function [e] = Set_2_e(e,m,Vc)
            % >> Error distribution.
            %  > e_c.
            e.c.c     = m.At\e.t.c;
            %  > abs().
            e.c.c_abs = abs(e.c.c);
            e.t.c_abs = abs(e.t.c);
            for i = 1:size(e.t.f,2)
                e.t.f_abs{i} = abs(e.t.f{i});
            end
            % >> Error norms.
        end 
    end
end
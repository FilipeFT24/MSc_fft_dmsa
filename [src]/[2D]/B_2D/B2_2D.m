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
                        k = ns;
                    otherwise
                        k = ns-1;
                end
                for j = 1:k
                    e.(i){j}.t.c       = zeros(Nc,1);
                    e.(i){j}.t.c_abs   = zeros(Nc,1);
                    e.(i){j}.t.f       = zeros(Nf,nc+1);
                    e.(i){j}.t.f_abs   = zeros(Nf,nc+1);
                    e.(i){j}.t.n.c     = zeros(1,3);
                    e.(i){j}.t.n_abs.c = zeros(1,3);
                    e.(i){j}.t.n.f     = zeros(3,nc+1);
                    e.(i){j}.t.n_abs.f = zeros(3,nc+1);
                    e.(i){j}.c.c       = zeros(Nc,1);
                    e.(i){j}.c.c_abs   = zeros(Nc,1);
                    e.(i){j}.c.n       = zeros(1,3);
                    e.(i){j}.c.n_abs   = zeros(1,3);
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update error.
        function [e] = Update_e(inp,msh,e,f,m,s,u,x)
            
           %e = B2_2D.Set_1_e(msh,m{1},e,f,inp.c,x{1}.xf.a,x{1}.xf.x,x{1}.nv.a.f,x{1}.nv.x.c);
            
            
        end
        % >> 1.3. ---------------------------------------------------------
        %  > 1.3.1. -------------------------------------------------------
        %  > Auxiliary function #1.
        %  > Compute remaining error fields (based on the convective/diffusive facial components).
        function [e] = Set_1_e(msh,m,e,f,v,x,z,y,xc)
            %  > \tau_f(\phi) and \tau(\nabla\phi).
            %etf
            
            
            etf(:,1) = v(1,1).*(x{1}-z{1})+v(2,1).*(x{2}(:,1)-z{2}(:,1)); % x
            etf(:,2) = v(1,2).*(x{1}-z{1})+v(2,2).*(x{2}(:,2)-z{2}(:,2)); % y
            
            etc = zeros(numel(xc),1);
            
            for i = 1:msh.c.Nc
                for j = 1:numel(msh.c.f.if(i,:))
                    k = msh.c.f.if(i,j);
                    etc_aux{i}(j,1) = msh.c.f.Sf{i}(j,:)*etf(k,:)';
                end
                etc(i,1) = sum(etc_aux{i});
            end
            
            
            l = f.av.c;
            e_1 = etc;
            e_2 = m.At*(l-xc);
          
            hold on;
            plot(e_1,'-or');
            plot(e_2,'-b');
            
            
        end
    end
end
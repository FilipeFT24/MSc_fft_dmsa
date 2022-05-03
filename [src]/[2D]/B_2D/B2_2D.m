classdef B2_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Update matrices.
        function [m] = Update_3(inp,msh,f,m,s,u,x)
            m        = B2_2D.Update_m(inp,msh,f,m,s,u,x);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update nodal/face values, etc.
        function [x] = Update_4(f,s,u,x)
            x        = Tools.Update_xv(f,s,u,x);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
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
    end
end
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
                        f_jk           = msh.c.f.if(j,k);
                        %  > Logical array w/ cell(s)/faces(s).
                        i_jk           = s.i      {f_jk,i};
                        l           = s.logical{f_jk,i};
                        %  > Cell/face indices used to fit "f_jk".
                        vA             = i_jk     ( l);
                        vB             = i_jk     (~l);
                        %  > Ac (cell contributions).
                        m.Ac {i}(j,vA) = m.Ac{i}(j,vA)+msh.c.f.Sf{j}(k,:)*x.Tf_V{f_jk,i}(:,l);
                        %  > Bc (cell contributions).
                        if any(~l)
                            fi           = ismembc(f.bd.i,sort(vB));
                            fv           = f.bd.v (fi);
                            m.Bc{i}(j,1) = m.Bc{i}(j)-msh.c.f.Sf{j}(k,:)*x.Tf_V{f_jk,i}(:,~l)*fv;
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
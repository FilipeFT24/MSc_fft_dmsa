classdef B_1D
    methods(Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [obj,pde] = Initialize(inp,msh)
            %  > Auxiliary variables.
            n  = length(inp.pv.v);
            m  = inp.pv.ns;
            Nc = msh.c.Nc;
            Nf = msh.f.Nf;
            
            % >> 'pde'.
            pde = B_1_1D.Update_pde(inp,msh);
            % >> 'obj'.
            %  > Field: 'e.a' (error w/ analytic field).
            obj.e.a.t.c       = zeros(Nc,1);
            obj.e.a.t.c_abs   = zeros(Nc,1);
            obj.e.a.t.f       = zeros(Nf,n+1);
            obj.e.a.t.f_abs   = zeros(Nf,n+1);
            obj.e.a.t.n.c     = zeros(1,3);
            obj.e.a.t.n_abs.c = zeros(1,3);
            obj.e.a.t.n.f     = zeros(3,n+1);
            obj.e.a.t.n_abs.f = zeros(3,n+1);
            obj.e.a.c.c       = zeros(Nc,1);
            obj.e.a.c.c_abs   = zeros(Nc,1);
            obj.e.a.c.n       = zeros(1,3);
            obj.e.a.c.n_abs   = zeros(1,3);  
            %  > Field: 'e.p' (error w/ predicted field).
            for i = 1:m
                obj.e.p{i}.t.c       = zeros(Nc,1);
                obj.e.p{i}.t.c_abs   = zeros(Nc,1);
                obj.e.p{i}.t.f       = zeros(Nf,n+1);
                obj.e.p{i}.t.f_abs   = zeros(Nf,n+1);
                obj.e.p{i}.t.n.c     = zeros(1,3);
                obj.e.p{i}.t.n_abs.c = zeros(1,3);
                obj.e.p{i}.t.n.f     = zeros(3,n+1);
                obj.e.p{i}.t.n_abs.f = zeros(3,n+1);
                obj.e.p{i}.c.c       = zeros(Nc,1);
                obj.e.p{i}.c.c_abs   = zeros(Nc,1);
                obj.e.p{i}.c.n       = zeros(1,3);
                obj.e.p{i}.c.n_abs   = zeros(1,3); 
            end
            %  > Field: 'm' (matrices).
            for i = 1:n
                obj.m.Ac{i} = zeros(Nc);
                obj.m.Af{i} = zeros(Nf,Nc);
                obj.m.Bc{i} = zeros(Nc,1);
                obj.m.Bf{i} = zeros(Nf,1);
            end
            obj.m.At        = zeros(Nc);
            obj.m.Bt        = pde.fn.vol;
            %  > Field: 's' (stencil coordinates,etc.).
            obj.s.c         = cell (Nf,n);
            obj.s.t         = cell (Nf,n);
            obj.s.bt        = cell (Nf,n);
            obj.s.bv        = cell (Nf,n);
            obj.s.v         = [inp.pv.v(1),-inp.pv.v(2)];
            %  > Field: 'u' (update stencil).
            obj.u           = A_2_1D.Initialize_upd(msh,inp.ps.p,inp.ps.t);
            %  > Field: 'x' (nodal solution/ stencil coefficients,etc.).
            obj.x.nv.a.c    = zeros(Nc,1);
            obj.x.nv.a.f    = zeros(Nf,1);
            obj.x.nv.x.c    = zeros(Nc,1);
            obj.x.cf        = cell (Nf,n);
            obj.x.vf.a      = cell (Nf,n);
            obj.x.vf.x      = cell (Nf,n);
            obj.x.xf.a      = zeros(Nf,n);
            obj.x.xf.x      = zeros(Nf,n);
            obj.x.if        = cell (Nf,n);                      
        end
        % >> 1.2. ---------------------------------------------------------
        function [obj,msh] = Run_p(inp,msh)
            %  > Initialize problem.
            [obj,pde] = B_1D.Initialize(inp,msh);
            
            %  > Set up 'p-standard' and 'p-adaptative' runs.
            switch inp.pa.adapt
                case false
                    %  > 'p-standard' run.
                    [obj,msh] = B_2_2_1D.p_standard(inp,obj,msh,pde);
                    %  > Plot...
                    if inp.pl.all
                        Fig_1_1D.Plot(obj,msh);
                    end
                case true
                    %  > 'p-standard' run.
                    [obj,msh] = B_2_2_1D.p_adaptative(inp,obj,msh,pde);
                    %  > Set nnz per cycle...
                    obj.e     = B_1D.Set_nnz(obj.e,obj.m); 
                otherwise
                    return;
            end
        end
        %  > 1.2.1. -------------------------------------------------------
        %  > Set nnz per cycle...
        function [e] = Set_nnz(e,m)
            for i = 1:size(m,1)
                for j = 1:size(m(i).Ac,2)
                    e{i,1}.nnz(j) = nnz(m(i).Ac{j});
                end
            end
        end
    end
end
classdef B_1D
    methods(Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [obj] = Initialize(inp,msh)
            %  > Auxiliary variables.
            n  = length(inp.pv.v);
            Nc = msh.c.Nc;
            Nf = msh.f.Nf;
            v  = inp.pv.v(1);
            g  = inp.pv.v(2);
                 
            %  > Fields: 'e.a', 'e.d' and 'e.p'.
            for i = ["a","da","p","d"]
                switch i
                    case "a"
                        ns = inp.pa.ns+1;
                    otherwise
                        ns = inp.pa.ns;
                end
                for j = 1:ns
                    obj.e.(i){j}.t.c       = zeros(Nc,1);
                    obj.e.(i){j}.t.c_abs   = zeros(Nc,1);
                    obj.e.(i){j}.t.f       = zeros(Nf,n+1);
                    obj.e.(i){j}.t.f_abs   = zeros(Nf,n+1);
                    obj.e.(i){j}.t.n.c     = zeros(1,3);
                    obj.e.(i){j}.t.n_abs.c = zeros(1,3);
                    obj.e.(i){j}.t.n.f     = zeros(3,n+1);
                    obj.e.(i){j}.t.n_abs.f = zeros(3,n+1);
                    obj.e.(i){j}.c.c       = zeros(Nc,1);
                    obj.e.(i){j}.c.c_abs   = zeros(Nc,1);
                    obj.e.(i){j}.c.n       = zeros(1,3);
                    obj.e.(i){j}.c.n_abs   = zeros(1,3);
                end
            end
            %  > Field: 'f' (analytic function).
            obj.f = B_1_1D.Update_func(inp,msh);
            %  > Field: 'm' (matrices).
            for i = 1:n
                obj.m.Ac{i} = zeros(Nc);
                obj.m.Af{i} = zeros(Nf,Nc);
                obj.m.Bc{i} = zeros(Nc,1);
                obj.m.Bf{i} = zeros(Nf,1);
            end
            obj.m.At     = zeros(Nc);
            obj.m.Bt     = obj.f.fv;
            %  > Field: 's' (stencil coordinates,etc.).
            obj.s.c      = cell (Nf,n);
            obj.s.t      = cell (Nf,n);
            obj.s.v      = [v,-g];
            obj.s.gv     = g./v;
            %  > Field: 'u' (update stencil).
            obj.u        = A_2_1D.Initialize_upd(msh,inp.pv.p);
            %  > Field: 'x' (nodal solution/ stencil coefficients,etc.).
            obj.x.if     = cell (Nf,n);
            obj.x.Tf     = cell (Nf,n);
            obj.x.nv.a.f = zeros(Nf,1);
            for i = ["a","x"]
                obj.x.cf.(i)   = cell (Nf,n);   %  > Unnecessary.
                obj.x.vf.(i)   = cell (Nf,n);
                obj.x.xf.(i)   = zeros(Nf,n);
                obj.x.nv.(i).c = zeros(Nc,1);
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Set up 'p-standard' and 'p-adaptative' runs.
        function [obj] = Run_p(inp,msh)
            switch inp.pa.adapt
                case false
                    %  > 'p-standard' run.
                    obj_p = B_1D.Initialize    (inp,msh);
                    obj   = B_2_2_1D.p_standard(inp,msh,obj_p);                    
                    %  > Plot...
                    Fig_V1_1_1D.Plot(inp,msh,obj,0);
                    Fig_V1_2_1D.Plot(inp,msh,obj,1);
                case true
                    %  > 'p-adaptative' run.
                    fld_adapt   = "p";
                    obj_adapt   = B_1D.Initialize(inp,msh);
                    for i       = 1:length(fld_adapt)
                        j       = fld_adapt(i);
                        obj.(j) = B_2_2_1D.p_adaptative(inp,obj_adapt,msh,j);
                    end
                    %  > Plot...
                    Fig_V1_3_1D.Plot(inp,obj,msh);
                otherwise
                    return;
            end
        end
    end
end
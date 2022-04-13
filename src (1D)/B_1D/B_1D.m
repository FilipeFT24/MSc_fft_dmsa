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
            f = ["a","p","d"];
            for i = 1:length(f)
                j = f(i);
                for k = 1:inp.pa.ns
                    obj.e.(j){k}.t.c       = zeros(Nc,1);
                    obj.e.(j){k}.t.c_abs   = zeros(Nc,1);
                    obj.e.(j){k}.t.f       = zeros(Nf,n+1);
                    obj.e.(j){k}.t.f_abs   = zeros(Nf,n+1);
                    obj.e.(j){k}.t.n.c     = zeros(1,3);
                    obj.e.(j){k}.t.n_abs.c = zeros(1,3);
                    obj.e.(j){k}.t.n.f     = zeros(3,n+1);
                    obj.e.(j){k}.t.n_abs.f = zeros(3,n+1);
                    obj.e.(j){k}.c.c       = zeros(Nc,1);
                    obj.e.(j){k}.c.c_abs   = zeros(Nc,1);
                    obj.e.(j){k}.c.n       = zeros(1,3);
                    obj.e.(j){k}.c.n_abs   = zeros(1,3);
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
            obj.x.Df     = cell (Nf,n);
            obj.x.if     = cell (Nf,n);
            obj.x.Tf     = cell (Nf,n);
            obj.x.nv.a.c = zeros(Nc,1);
            obj.x.nv.a.f = zeros(Nf,1);
            obj.x.nv.x.c = zeros(Nc,1);
            obj.x.cf.a   = cell (Nf,n);
            obj.x.cf.x   = cell (Nf,n);
            obj.x.vf.a   = cell (Nf,n);
            obj.x.vf.x   = cell (Nf,n);
            obj.x.xf.a   = zeros(Nf,n);   %  > Face  nodal values (w/ '.a').
            obj.x.xf.x   = zeros(Nf,n);   %  > Face  nodal values (w/ '.x').
            obj.x.xt.a   = cell (Nf,n);   %  > Total nodal values (w/ '.a').
            obj.x.xt.x   = cell (Nf,n);   %  > Total nodal values (w/ '.x').
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
                    plot  = [1,0];
                    if plot(1)
                        Fig_V1_1_1D.Plot(inp,msh,obj);
                    end
                    if plot(2)
                        Fig_V1_2_1D.Plot(inp,msh,obj);
                    end
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
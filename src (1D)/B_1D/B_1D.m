classdef B_1D
    methods(Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [obj,pde] = Initialize(inp,msh)
            %  > Auxiliary variables.
            n  = length(inp.pv.v);
            Nc = msh.c.Nc;
            Nf = msh.f.Nf;
            
            % >> 'pde'.
            pde = B_1_1D.Update_pde(inp,msh);
            % >> 'obj'.
            %  > Fields: 'e.a', 'e.d', 'e.p' and 'e.t'.
            f = ["a","p","d","t"];
            for i = 1:length(f)
                %  > Fields: 'a' and 'p'.
                if i <= 2
                    if i == 1
                        if ~inp.pa.comp_av
                            m = 2;
                        else
                            m = inp.pa.ns+1;
                        end
                    else
                        m = inp.pa.ns;
                    end
                    for j = 1:m
                        obj.e.(f(i)){j}.t.c       = zeros(Nc,1);
                        obj.e.(f(i)){j}.t.c_abs   = zeros(Nc,1);
                        obj.e.(f(i)){j}.t.f       = zeros(Nf,n+1);
                        obj.e.(f(i)){j}.t.f_abs   = zeros(Nf,n+1);
                        obj.e.(f(i)){j}.t.n.c     = zeros(1,3);
                        obj.e.(f(i)){j}.t.n_abs.c = zeros(1,3);
                        obj.e.(f(i)){j}.t.n.f     = zeros(3,n+1);
                        obj.e.(f(i)){j}.t.n_abs.f = zeros(3,n+1);
                        obj.e.(f(i)){j}.c.c       = zeros(Nc,1);
                        obj.e.(f(i)){j}.c.c_abs   = zeros(Nc,1);
                        obj.e.(f(i)){j}.c.n       = zeros(1,3);
                        obj.e.(f(i)){j}.c.n_abs   = zeros(1,3);
                    end
                end
                %  > Field: 'd'.
                if i == 3
                    for j = 1:inp.pa.ns
                        obj.e.(f(i)){j}.t.f       = zeros(Nf,n);
                        obj.e.(f(i)){j}.t.f_abs   = zeros(Nf,n);
                        obj.e.(f(i)){j}.t.n.f     = zeros(3,n);
                        obj.e.(f(i)){j}.t.n_abs.f = zeros(3,n);
                    end
                end
                %  > Field 't'.
                if i == 4 && inp.pt.tt
                    obj.e.(f(i)).c       = cell(1,n);
                    obj.e.(f(i)).c_abs   = cell(1,n);
                    obj.e.(f(i)).f       = cell(1,n);
                    obj.e.(f(i)).f_abs   = cell(1,n);
                end
            end
            %  > Field: 'm' (matrices).
            for i = 1:n
                obj.m.Ac{i} = zeros(Nc);
                obj.m.Af{i} = zeros(Nf,Nc);
                obj.m.Bc{i} = zeros(Nc,1);
                obj.m.Bf{i} = zeros(Nf,1);
            end
            obj.m.At     = zeros(Nc);
            obj.m.Bt     = pde.fn.vol;
            %  > Field: 's' (stencil coordinates,etc.).
            obj.s.c      = cell (Nf,n);
            obj.s.t      = cell (Nf,n);
            obj.s.bt     = cell (Nf,n);
            obj.s.bv     = cell (Nf,n);
            obj.s.v      = [inp.pv.v(1),-inp.pv.v(2)];
            %  > Field: 'u' (update stencil).
            obj.u        = A_2_1D.Initialize_upd(msh,inp.ps.p,inp.ps.t);
            %  > Field: 'x' (nodal solution/ stencil coefficients,etc.).
            obj.x.nv.a.c = zeros(Nc,1);
            obj.x.nv.a.f = zeros(Nf,1);
            obj.x.nv.x.c = zeros(Nc,1);
            obj.x.cf     = cell (Nf,n);
            obj.x.vf.a   = cell (Nf,n);
            obj.x.vf.x   = cell (Nf,n);
            obj.x.xf.a   = zeros(Nf,n);
            obj.x.xf.x   = zeros(Nf,n);
            obj.x.if     = cell (Nf,n);
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
                    Fig_V1_1_1D.Plot(obj,msh);
                    if inp.pt.tt
                        Fig_V1_2_1D.Plot(inp,obj,msh);
                    end
                case true
                    %  > 'p-standard' run.
                    [obj,msh] = B_2_2_1D.p_adaptative(inp,obj,msh,pde);
                    %  > Plot...
                    Fig_V1_1_1D.Plot(obj,msh);
                otherwise
                    return;
            end
        end
    end
end
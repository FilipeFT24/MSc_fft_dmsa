classdef Tools_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [] = Set_Directories()
            addpath(genpath('A_1D'));
            addpath(genpath('B_1D'));
            addpath(genpath('C_1D'));
            addpath(genpath('../[Post-processing]'));
            addpath(genpath('../[Tools]')); 
        end
        % >> 1.2. ---------------------------------------------------------
        function [obj]  = Sort_obj(obj)
            obj         = orderfields(obj        ,{'e','m','s','u','x'});
            obj.e       = orderfields(obj.e      ,{'a','p'});
            obj.e.a     = orderfields(obj.e.a    ,{'c','t'});
            obj.e.a.c   = orderfields(obj.e.a.c  ,{'c','c_abs','n','n_abs'});
            obj.e.a.t   = orderfields(obj.e.a.t  ,{'c','c_abs','f','f_abs','n','n_abs'});
            obj.e.a.t.n = orderfields(obj.e.a.t.n,{'c','f'});
            obj.m       = orderfields(obj.m      ,{'Ac','Af','At','Bc','Bf','Bt','nnz'});
            obj.s       = orderfields(obj.s      ,{'bt','bv','c','t','v'});
            obj.u       = orderfields(obj.u      ,{'p','s'});
            obj.x       = orderfields(obj.x      ,{'cf','if','nv','vf','xf'});
        end
        % >> 1.3. ---------------------------------------------------------
        function [msh] = Sort_msh(msh)
            msh     = orderfields(msh    ,{'c','d','f'});
            msh.c   = orderfields(msh.c  ,{'f','Nc','Vc','Xc'});
            msh.c.f = orderfields(msh.c.f,{'f','Nf'});
            msh.f   = orderfields(msh.f  ,{'c','Nf','Xv'});
        end
        % >> 1.4. ---------------------------------------------------------
        function [pde] = Sort_pde(pde)
            pde     = orderfields(pde   ,{'av','fn'});
            pde.av  = orderfields(pde.av,{'c','f'});
            pde.fn  = orderfields(pde.fn,{'f','i','vol'});
        end

        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Compute error norms (cell/face L1,L2 and L_infinity norms).
        function [L] = Set_n(E,V)
            if nargin == 1
                L(1,:) = mean(E);
                L(2,:) = mean(sqrt(E.^2));
                L(3,:) = max(E);
            else
                L(1,:) = sum (E.*V)./sum(V);
                L(2,:) = sum (sqrt((E.*V).^2))./sum(sqrt(V.^2));
                L(3,:) = max(E);
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Compute remaining error fields (based on the convective/diffusive facial components).
        function [e] = Set_e(e,m,Vc)
            % >> Error distribution.
            %  > \tau_f.
            [a,b]          = size(e.t.f);
            e.t.f    (:,b) = sum (e.t.f(:,1:b-1),2);
            %  > \tau_c.
            for i = 1:a-1
                e.t.c(i,1) = e.t.f(i,b)-e.t.f(i+1,b);
            end
            %  > e_c.
            e.c.c    (:,1) = m.At\e.t.c;
            %  > abs().
            e.c.c_abs      = abs(e.c.c);
            e.t.c_abs      = abs(e.t.c);
            e.t.f_abs      = abs(e.t.f);
            % >> Error norms.
            e.c.n          = Tools_1D.Set_n(e.c.c,Vc);
            e.c.n_abs      = Tools_1D.Set_n(e.c.c_abs,Vc);
            e.t.n.f        = Tools_1D.Set_n(e.t.f);
            e.t.n_abs.f    = Tools_1D.Set_n(e.t.f_abs);
            e.t.n.c        = Tools_1D.Set_n(e.t.c,Vc);
            e.t.n_abs.c    = Tools_1D.Set_n(e.t.c_abs,Vc);
        end
        % >> 2.3. ---------------------------------------------------------
        %  >
        function [] = pol_shp()
        end
    end
end
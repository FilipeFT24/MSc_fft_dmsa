classdef Tools_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [] = Set_Directories()
            addpath(genpath('A_1D'));
            addpath(genpath('B_1D'));
            addpath(genpath('C_1D'));
            addpath(genpath('D_1D'));
            addpath(genpath('../[Tools]/[Tools - Data]'));
            addpath(genpath('../[Tools]/[Tools - Mesh generation]'));
            addpath(genpath('../[Tools]/[Tools - Numerical]'));
            addpath(genpath('../[Tools]/[Tools - Other stuff]'));
        end
        % >> 1.2. ---------------------------------------------------------
        function [obj]  = Sort_obj(obj)
%             obj         = orderfields(obj        ,{'e','m','s','u','x'});
%             obj.e       = orderfields(obj.e      ,{'a','p'});
%             obj.e.a     = orderfields(obj.e.a    ,{'c','t'});
%             obj.e.a.c   = orderfields(obj.e.a.c  ,{'c','c_abs','n','n_abs'});
%             obj.e.a.t   = orderfields(obj.e.a.t  ,{'c','c_abs','f','f_abs','n','n_abs'});
%             obj.e.a.t.n = orderfields(obj.e.a.t.n,{'c','f'});
%             obj.m       = orderfields(obj.m      ,{'Ac','Af','At','Bc','Bf','Bt','nnz'});
%             obj.s       = orderfields(obj.s      ,{'bt','bv','c','t','v'});
%             obj.u       = orderfields(obj.u      ,{'p','s'});
%             obj.x       = orderfields(obj.x      ,{'cf','if','nv','vf','xf'});
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
        function [L] = n(E,V)
            if nargin == 1
                L(1,:) = Tools_1D.L1(E);
                L(2,:) = Tools_1D.L2(E);
            else
                L(1,:) = Tools_1D.L1(E,V);
                L(2,:) = Tools_1D.L2(E,V);
            end
            L(3,:) = Tools_1D.L3(E);
        end
        %  > L1.
        function [L1] = L1(E,V)
            if nargin == 1
                L1 = mean(E);
            else
                L1 = sum (E.*V)./sum(V);
            end
        end
        %  > L2.
        function [L2] = L2(E,V)
            if nargin == 1
                L2 = mean(sqrt(E.^2));
            else
                L2 = sum (sqrt((E.*V).^2))./sum(sqrt(V.^2));
            end
        end
        %  > L_infinity.
        function [L3] = L3(E)
            L3 = max(E);
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Matrix inversion w/ 'p-adaptative' process (...to compute "updated" cell global discretization error).
        function [y] = p_adapt_inv(opt,str)
            switch opt
                case 1
                    a = str.a;
                    y = inv(a);
                case 2
                    a = str.a;
                    b = str.b;
                    i = eye(size(a));
                    y = a*(i-inv(i+b*a)*b*a);
                otherwise
                    return;
            end                   
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
    end
end
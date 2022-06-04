classdef Tools_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        
        % >> 1.2. ---------------------------------------------------------
        %  > 1.2.1. -------------------------------------------------------
        function [msh] = Sort_msh(msh)
            msh     = orderfields(msh    ,{'c','d','f'});
            msh.c   = orderfields(msh.c  ,{'f','Nc','Vc','Xc'});
            msh.c.f = orderfields(msh.c.f,{'f','Nf'});
            msh.f   = orderfields(msh.f  ,{'c','Nf','Xv'});
        end
        %  > 1.2.2. -------------------------------------------------------
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
    end
end
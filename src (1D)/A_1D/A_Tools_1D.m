classdef A_Tools_1D 
    methods (Static)
        %% > A_Tools (1D).
        % >> --------------------------------------------------------------
        % >> 1.   Sort 'msh' fields.
        % >> --------------------------------------------------------------
        
        %% > 1. -----------------------------------------------------------      
        % >> 1.2. ---------------------------------------------------------
        function [msh] = Sort_msh(msh)
            % >> msh.
            msh     = orderfields(msh        ,{'d','c','f','bnd','s'});
            %  > c.
            msh.c   = orderfields(msh.c      ,{'NC','xy_v','mean','h','vol','v','c','f'});
            msh.c.f = orderfields(msh.c.f    ,{'f','xy_v','mean','len','Nf','Sf'});
            %  > f.
            msh.f   = orderfields(msh.f      ,{'NF','xy_v','mean','v','c'});
            %  > bnd.
            msh.bnd = orderfields(msh.bnd    ,{'c','f'});
            %  > s.
            msh.s     = orderfields(msh.s    ,{'c','f','c_e','f_e','xy_v_c','xy_v_f','xy_v_t','par'});
            msh.s.par = orderfields(msh.s.par,{'ne','nx','ny','lx','ly'});
        end
    end
end
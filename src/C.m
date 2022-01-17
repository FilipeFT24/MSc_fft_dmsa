classdef C
    methods(Static)
        %% > Wrap-up C.
        % >> --------------------------------------------------------------
        % >> 1. Figure 0.
        % >> 2. Figure 1.
        % >> --------------------------------------------------------------
        function [] = WrapUp_C(inp,msh,blk,dt,len)
            %  > Figure 0.
            %  Fig_0.WrapUp_Fig_0(1,inp.fr.ng);
            %  > Figure 1.
            Fig_1.WrapUp_Fig_1(2,inp,msh,dt,len);
            %  > Figure 2.
            Fig_2.WrapUp_Fig_2(3,inp,msh,blk);
        end
    end
end
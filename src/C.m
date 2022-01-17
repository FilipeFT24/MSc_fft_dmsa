classdef C
    methods(Static)
        %% > Wrap-up C.
        % >> --------------------------------------------------------------
        % >> 1. Figure 0.
        % >> 2. Figure 1.
        % >> --------------------------------------------------------------
        function [] = WrapUp_C(inp,msh,pde,dt,len)
            %  > Working directory.
            Tools.Set_Directory('C');
            %  > Figure 0.
            %  Fig_0.WrapUp_Fig_0(1,inp.fr.ng);
            %  > Figure 1.
            Fig_1.WrapUp_Fig_1(2,inp,msh,pde,dt,len);
            %  > Figure 2.
            %  Fig_2.WrapUp_Fig_2(3,inp,msh,pde);
        end
    end
end
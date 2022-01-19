classdef C
    methods(Static)
        %% > Wrap-up C.
        function [] = WrapUp_C(inp,msh,pde,dt,len)
            %  > Working directory.
            Tools.Set_Directory('C');

            %  > Figure 0.
            Fig_0.WrapUp_Fig_0(false,false,1,inp.fr.ng);
            %  > Figure 1.
            Fig_1.WrapUp_Fig_1(true ,true,2,inp,msh,dt,len);
            %  > Figure 2.
            Fig_2.WrapUp_Fig_2(false,3,inp,msh,pde,len);
            %  > Figure 3.
            %  > Figure 4.
            Fig_4.WrapUp_Fig_4(false,5,inp,msh,pde,len);
        end
    end
end
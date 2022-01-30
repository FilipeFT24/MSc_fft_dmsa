classdef C
    methods(Static)
        %% > Wrap-up C.
        function [] = WrapUp_C(inp,msh,pde,dt,len)
            %  > Figure 0.
            Fig_0.WrapUp_Fig_0(false,false,1,inp.fr.ng);
            %  > Figure 1.
            Fig_1.WrapUp_Fig_1(true ,false,2,inp,msh,dt,len);
            %  > Figure 2.
            Fig_2.WrapUp_Fig_2(false,false,3,inp,msh,pde,len)
        end
    end
end
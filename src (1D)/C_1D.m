classdef C_1D
    methods(Static)
        %% > Wrap-up C (1D).
        function [] = WrapUp_C_1D(msh,pde)
            %  > Figure 1.
            Fig_1_1D.WrapUp_Fig_1_1D(false,false,1,msh);
            %  > Figure 2.
            Fig_2_1D.WrapUp_Fig_2_1D(true ,false,2,msh,pde);
        end
    end
end